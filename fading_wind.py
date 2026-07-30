#!/usr/bin/env python3
"""Estimate horizontal wind from three-dipole power-fading correlations."""

from __future__ import annotations

import numpy as np
import scipy.optimize as optimize
import scipy.signal as signal


CHANNEL_FIELDS = ("rti1", "rti3", "rti4")
# WGS84 east/north coordinates relative to channel 3 (MF1).
ANTENNA_EN_M = np.asarray(
    (
        (144.5, -79.2),  # ch1, MF3
        (0.0, 0.0),      # ch3, MF1
        (4.7, -165.1),   # ch4, MF2
    ),
    dtype=np.float64,
)
PAIR_INDICES = ((0, 1), (0, 2), (1, 2))
ALTITUDE_KM = np.arange(70.0, 125.0, 5.0)
ALTITUDE_HALF_WIDTH_KM = 2.0
WINDOW_S = 5 * 60
OUTPUT_CADENCE_S = 60
MAX_LAG_S = 30.0
MIN_PEAK_CORRELATION = 0.15
MAX_FIT_RMSE = 0.35
MAX_PATTERN_SPEED_MS = 600.0
CORRELATION_SAMPLE_RATE_HZ = 4.0


def _decimate_power(
    power: np.ndarray,
    sample_interval_s: float,
) -> tuple[np.ndarray, float]:
    decimation = max(
        1,
        int(round(1.0 / (CORRELATION_SAMPLE_RATE_HZ * sample_interval_s))),
    )
    if decimation == 1:
        return power, sample_interval_s
    usable = len(power) // decimation * decimation
    return (
        power[:usable].reshape(
            usable // decimation,
            decimation,
            power.shape[1],
        ).mean(axis=1),
        sample_interval_s * decimation,
    )


def _normalized_correlation(
    first: np.ndarray,
    second: np.ndarray,
    maximum_lag_samples: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return corr(first(t), second(t+lag)) and integer sample lags."""
    first = signal.detrend(np.asarray(first, dtype=np.float64), type="linear")
    second = signal.detrend(np.asarray(second, dtype=np.float64), type="linear")
    first_scale = np.std(first)
    second_scale = np.std(second)
    if first_scale <= 0.0 or second_scale <= 0.0:
        raise ValueError("constant fading series")
    first /= first_scale
    second /= second_scale
    correlation = signal.correlate(second, first, mode="full", method="fft")
    lags = signal.correlation_lags(len(second), len(first), mode="full")
    keep = np.abs(lags) <= maximum_lag_samples
    lags = lags[keep]
    overlap = len(first) - np.abs(lags)
    return correlation[keep] / np.maximum(overlap, 1), lags


def _initial_pattern_velocity(
    correlations: list[np.ndarray],
    lag_seconds: np.ndarray,
) -> np.ndarray:
    peak_lags = np.asarray(
        [lag_seconds[int(np.nanargmax(value))] for value in correlations]
    )
    baselines = np.asarray(
        [
            ANTENNA_EN_M[second] - ANTENNA_EN_M[first]
            for first, second in PAIR_INDICES
        ]
    )
    inverse_velocity, *_ = np.linalg.lstsq(
        baselines,
        peak_lags,
        rcond=None,
    )
    denominator = float(inverse_velocity @ inverse_velocity)
    if denominator <= 1e-12:
        return np.zeros(2, dtype=np.float64)
    velocity = inverse_velocity / denominator
    speed = np.linalg.norm(velocity)
    if speed > MAX_PATTERN_SPEED_MS:
        velocity *= MAX_PATTERN_SPEED_MS / speed
    return velocity


def cross_correlation_delays(
    power: np.ndarray,
    sample_interval_s: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return sub-sample peak delay and peak coefficient for each baseline."""
    power, sample_interval_s = _decimate_power(
        np.asarray(power, dtype=np.float64),
        sample_interval_s,
    )
    maximum_lag_samples = int(round(MAX_LAG_S / sample_interval_s))
    delays_s = np.full(len(PAIR_INDICES), np.nan)
    peak_correlation = np.full(len(PAIR_INDICES), np.nan)
    for pair_index, (first, second) in enumerate(PAIR_INDICES):
        correlation, lag_samples = _normalized_correlation(
            power[:, first],
            power[:, second],
            maximum_lag_samples,
        )
        peak = int(np.nanargmax(correlation))
        fractional_peak = 0.0
        if 0 < peak < len(correlation) - 1:
            lower = correlation[peak - 1]
            center = correlation[peak]
            upper = correlation[peak + 1]
            denominator = lower - 2.0 * center + upper
            if abs(denominator) > 1e-12:
                fractional_peak = float(
                    np.clip(
                        0.5 * (lower - upper) / denominator,
                        -1.0,
                        1.0,
                    )
                )
        delays_s[pair_index] = (
            lag_samples[peak] + fractional_peak
        ) * sample_interval_s
        peak_correlation[pair_index] = correlation[peak]
    return delays_s, peak_correlation


def fit_fading_window(
    power: np.ndarray,
    sample_interval_s: float,
) -> tuple[float, float, float, float]:
    """
    Fit an evolving elliptical drifting-pattern correlation model.

    Returns zonal wind, meridional wind, median cross-correlation peak, and
    normalized correlation-fit RMSE. The diffraction-pattern velocity is
    divided by two to obtain neutral wind.
    """
    power = np.asarray(power, dtype=np.float64)
    if power.ndim != 2 or power.shape[1] != len(CHANNEL_FIELDS):
        raise ValueError("power must have shape (time, three dipoles)")
    power, sample_interval_s = _decimate_power(power, sample_interval_s)
    maximum_lag_samples = int(round(MAX_LAG_S / sample_interval_s))
    if len(power) < 4 * maximum_lag_samples:
        raise ValueError("fading window is too short")

    cross_correlations = []
    lag_samples = None
    for first, second in PAIR_INDICES:
        correlation, current_lags = _normalized_correlation(
            power[:, first],
            power[:, second],
            maximum_lag_samples,
        )
        cross_correlations.append(correlation)
        if lag_samples is None:
            lag_samples = current_lags
    auto_correlations = [
        _normalized_correlation(
            power[:, channel],
            power[:, channel],
            maximum_lag_samples,
        )[0]
        for channel in range(len(CHANNEL_FIELDS))
    ]
    lag_seconds = lag_samples.astype(np.float64) * sample_interval_s
    observed = np.asarray(
        [np.mean(auto_correlations, axis=0), *cross_correlations]
    )
    baselines = np.asarray(
        [
            (0.0, 0.0),
            *(
                ANTENNA_EN_M[second] - ANTENNA_EN_M[first]
                for first, second in PAIR_INDICES
            ),
        ],
        dtype=np.float64,
    )
    initial_velocity = _initial_pattern_velocity(
        cross_correlations,
        lag_seconds,
    )

    def model_and_residual(parameters: np.ndarray) -> np.ndarray:
        velocity = parameters[:2]
        cholesky = np.asarray(
            (
                (np.exp(parameters[2]), 0.0),
                (parameters[3], np.exp(parameters[4])),
            )
        ) / 200.0
        precision = cholesky @ cholesky.T
        inverse_decay_s = np.exp(parameters[5]) / 60.0
        model_rows = []
        for baseline, correlation in zip(baselines, observed):
            displacement = (
                baseline[None, :]
                - lag_seconds[:, None] * velocity[None, :]
            )
            exponent = -np.einsum(
                "ni,ij,nj->n",
                displacement,
                precision,
                displacement,
            ) - (lag_seconds * inverse_decay_s) ** 2
            shape = np.exp(np.clip(exponent, -50.0, 0.0))
            amplitude = np.clip(
                float(correlation @ shape)
                / max(float(shape @ shape), 1e-12),
                0.0,
                1.2,
            )
            model_rows.append(amplitude * shape)
        model = np.asarray(model_rows)
        return (model - observed).ravel()

    fit = optimize.least_squares(
        model_and_residual,
        np.asarray(
            (
                initial_velocity[0],
                initial_velocity[1],
                0.0,
                0.0,
                0.0,
                0.0,
            )
        ),
        bounds=(
            (-MAX_PATTERN_SPEED_MS, -MAX_PATTERN_SPEED_MS, -2.3, -5.0, -2.3, -2.3),
            (MAX_PATTERN_SPEED_MS, MAX_PATTERN_SPEED_MS, 2.3, 5.0, 2.3, 3.0),
        ),
        loss="soft_l1",
        f_scale=0.08,
        max_nfev=300,
    )
    pattern_velocity = fit.x[:2]
    pattern_speed = float(np.linalg.norm(pattern_velocity))
    peak_correlation = float(
        np.median([np.nanmax(value) for value in cross_correlations])
    )
    peak_delays_s = np.asarray(
        [
            lag_seconds[int(np.nanargmax(value))]
            for value in cross_correlations
        ]
    )
    fit_rmse = float(
        np.sqrt(np.mean(model_and_residual(fit.x) ** 2))
    )
    valid = (
        fit.success
        and peak_correlation >= MIN_PEAK_CORRELATION
        and fit_rmse <= MAX_FIT_RMSE
        and pattern_speed <= MAX_PATTERN_SPEED_MS
        and np.all(np.abs(peak_delays_s) < 0.95 * MAX_LAG_S)
    )
    if not valid:
        return np.nan, np.nan, peak_correlation, fit_rmse
    return (
        0.5 * float(pattern_velocity[0]),
        0.5 * float(pattern_velocity[1]),
        peak_correlation,
        fit_rmse,
    )


def load_altitude_power(
    reader,
    start_unix: float,
    end_unix: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Read 20 Hz metadata and average dipole power into altitude bands."""
    rows = []
    start_us = int(start_unix * 1e6)
    end_us = int(end_unix * 1e6)
    for chunk_start in range(start_us, end_us, int(60e6)):
        metadata = reader.read(
            chunk_start,
            min(chunk_start + int(60e6), end_us) - 1,
        )
        for key in sorted(metadata):
            record = metadata[key]
            required = (*CHANNEL_FIELDS, "rvec", "tvec")
            if not all(field in record for field in required):
                continue
            range_km = np.asarray(record["rvec"], dtype=np.float64)
            band_power = np.full(
                (
                    len(record["tvec"]),
                    len(CHANNEL_FIELDS),
                    len(ALTITUDE_KM),
                ),
                np.nan,
                dtype=np.float32,
            )
            for altitude_index, altitude_km in enumerate(ALTITUDE_KM):
                mask = (
                    (range_km >= altitude_km - ALTITUDE_HALF_WIDTH_KM)
                    & (range_km <= altitude_km + ALTITUDE_HALF_WIDTH_KM)
                )
                if not np.any(mask):
                    continue
                for channel_index, field in enumerate(CHANNEL_FIELDS):
                    band_power[:, channel_index, altitude_index] = np.mean(
                        np.abs(np.asarray(record[field])[:, mask]) ** 2,
                        axis=1,
                    )
            sample_times = (
                float(key) / 1e6
                + np.asarray(record["tvec"], dtype=np.float64)
            )
            rows.append((sample_times, band_power))
    if not rows:
        raise RuntimeError("no fading metadata in requested interval")
    return (
        np.concatenate([row[0] for row in rows]),
        np.concatenate([row[1] for row in rows], axis=0),
    )


def estimate_winds(
    reader,
    start_unix: float,
    end_unix: float,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Return minute-cadence zonal/meridional fading winds by altitude."""
    half_window = WINDOW_S / 2.0
    sample_times, power = load_altitude_power(
        reader,
        start_unix - half_window,
        end_unix,
    )
    order = np.argsort(sample_times)
    sample_times = sample_times[order]
    power = power[order]
    sample_interval_s = float(np.median(np.diff(sample_times)))
    first_center = int(
        np.ceil((start_unix + half_window) / OUTPUT_CADENCE_S)
        * OUTPUT_CADENCE_S
    )
    final_center = int(
        np.floor((end_unix - half_window) / OUTPUT_CADENCE_S)
        * OUTPUT_CADENCE_S
    )
    centers = np.arange(
        first_center,
        final_center + 1,
        OUTPUT_CADENCE_S,
        dtype=np.int64,
    )
    zonal = np.full((len(centers), len(ALTITUDE_KM)), np.nan)
    meridional = np.full_like(zonal, np.nan)
    peak_correlation = np.full_like(zonal, np.nan)
    fit_rmse = np.full_like(zonal, np.nan)
    baseline_delay_s = np.full(
        (len(centers), len(ALTITUDE_KM), len(PAIR_INDICES)),
        np.nan,
    )
    baseline_peak_correlation = np.full_like(baseline_delay_s, np.nan)
    for time_index, center in enumerate(centers):
        mask = (
            (sample_times >= center - half_window)
            & (sample_times < center + half_window)
        )
        expected = int(round(WINDOW_S / sample_interval_s))
        if np.count_nonzero(mask) < 0.9 * expected:
            continue
        for altitude_index in range(len(ALTITUDE_KM)):
            values = power[mask, :, altitude_index]
            if not np.all(np.isfinite(values)):
                continue
            (
                baseline_delay_s[time_index, altitude_index],
                baseline_peak_correlation[time_index, altitude_index],
            ) = cross_correlation_delays(values, sample_interval_s)
            (
                zonal[time_index, altitude_index],
                meridional[time_index, altitude_index],
                peak_correlation[time_index, altitude_index],
                fit_rmse[time_index, altitude_index],
            ) = fit_fading_window(values, sample_interval_s)
    return (
        centers,
        zonal,
        meridional,
        peak_correlation,
        fit_rmse,
        baseline_delay_s,
        baseline_peak_correlation,
    )
