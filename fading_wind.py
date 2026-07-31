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
ALTITUDE_KM = np.arange(50.0, 125.0, 5.0)
ALTITUDE_HALF_WIDTH_KM = 2.0
WINDOW_S = 5 * 60
OUTPUT_CADENCE_S = 60
MAX_LAG_S = 30.0
MIN_BASELINE_CORRELATION = 0.30
MIN_PEAK_PROMINENCE = 0.08
MIN_PEAK_UNIQUENESS = 0.04
MAX_FIT_RMSE = 0.18
MAX_MODEL_PEAK_ERROR_S = 1.0
MAX_NEUTRAL_WIND_SPEED_MS = 120.0
MAX_PATTERN_SPEED_MS = 2.0 * MAX_NEUTRAL_WIND_SPEED_MS
CORRELATION_SAMPLE_RATE_HZ = 4.0
TEMPORAL_BANDPASS_HZ = (1.0 / 120.0, 1.5)
SUBWINDOW_COUNT = 3
MIN_VALID_SUBWINDOWS = 2
MAX_SUBWINDOW_COMPONENT_SPREAD_MS = 50.0
MAX_BOOTSTRAP_UNCERTAINTY_MS = 25.0
BOOTSTRAP_REPLICATES = 500
MAX_COMMON_MODE_FRACTION = 0.85
MAX_SCALED_FIT_CONDITION = 50.0


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


def _condition_power(
    power: np.ndarray,
    sample_interval_s: float,
) -> np.ndarray:
    """Suppress slow common gain changes and impulsive meteor-like outliers."""
    power = np.asarray(power, dtype=np.float64)
    floor = np.maximum(
        np.nanpercentile(power, 1.0, axis=0),
        np.finfo(np.float64).tiny,
    )
    conditioned = np.log(np.maximum(power, floor[None, :]))
    median = np.nanmedian(conditioned, axis=0)
    mad = np.nanmedian(np.abs(conditioned - median[None, :]), axis=0)
    scale = np.maximum(1.4826 * mad, 1e-6)
    conditioned = np.clip(
        conditioned,
        median[None, :] - 8.0 * scale[None, :],
        median[None, :] + 8.0 * scale[None, :],
    )
    sample_rate_hz = 1.0 / sample_interval_s
    low_hz, high_hz = TEMPORAL_BANDPASS_HZ
    high_hz = min(high_hz, 0.45 * sample_rate_hz)
    if low_hz >= high_hz:
        raise ValueError("correlation sample rate is too low for band-pass")
    sos = signal.butter(
        3,
        (low_hz, high_hz),
        btype="bandpass",
        fs=sample_rate_hz,
        output="sos",
    )
    return signal.sosfiltfilt(sos, conditioned, axis=0)


def _peak_metrics(
    correlation: np.ndarray,
    lag_seconds: np.ndarray,
) -> dict[str, float]:
    peak_index = int(np.nanargmax(correlation))
    fractional_peak = 0.0
    if 0 < peak_index < len(correlation) - 1:
        lower = correlation[peak_index - 1]
        center = correlation[peak_index]
        upper = correlation[peak_index + 1]
        denominator = lower - 2.0 * center + upper
        if abs(denominator) > 1e-12:
            fractional_peak = float(
                np.clip(
                    0.5 * (lower - upper) / denominator,
                    -1.0,
                    1.0,
                )
            )
    lag_step_s = float(np.median(np.diff(lag_seconds)))
    delay_s = float(
        lag_seconds[peak_index] + fractional_peak * lag_step_s
    )
    if 0 < peak_index < len(correlation) - 1:
        prominence = float(
            signal.peak_prominences(correlation, [peak_index])[0][0]
        )
        width_samples = float(
            signal.peak_widths(
                correlation,
                [peak_index],
                rel_height=0.5,
            )[0][0]
        )
    else:
        prominence = 0.0
        width_samples = 0.0
    exclusion = max(2, int(np.ceil(max(width_samples, 1.0))))
    keep = np.ones(len(correlation), dtype=bool)
    keep[
        max(0, peak_index - exclusion):
        min(len(correlation), peak_index + exclusion + 1)
    ] = False
    secondary = (
        float(np.nanmax(correlation[keep]))
        if np.any(keep)
        else float(correlation[peak_index])
    )
    return {
        "delay_s": delay_s,
        "peak": float(correlation[peak_index]),
        "prominence": prominence,
        "width_s": width_samples * abs(lag_step_s),
        "uniqueness": float(correlation[peak_index] - secondary),
    }


def correlation_diagnostics(
    power: np.ndarray,
    sample_interval_s: float,
    *,
    conditioned: bool = True,
) -> dict[str, np.ndarray]:
    """Return full ACF/CCF curves and objective peak diagnostics."""
    power, sample_interval_s = _decimate_power(
        np.asarray(power, dtype=np.float64),
        sample_interval_s,
    )
    if conditioned:
        power = _condition_power(power, sample_interval_s)
    maximum_lag_samples = int(round(MAX_LAG_S / sample_interval_s))
    lag_samples = np.arange(
        -maximum_lag_samples,
        maximum_lag_samples + 1,
    )
    lag_seconds = lag_samples.astype(np.float64) * sample_interval_s
    cross_correlations = []
    cross_metrics = []
    for first, second in PAIR_INDICES:
        correlation, current_lags = _normalized_correlation(
            power[:, first],
            power[:, second],
            maximum_lag_samples,
        )
        if not np.array_equal(current_lags, lag_samples):
            raise RuntimeError("inconsistent correlation lag grid")
        cross_correlations.append(correlation)
        cross_metrics.append(_peak_metrics(correlation, lag_seconds))
    auto_correlations = [
        _normalized_correlation(
            power[:, channel],
            power[:, channel],
            maximum_lag_samples,
        )[0]
        for channel in range(len(CHANNEL_FIELDS))
    ]
    standardized = (
        power - np.mean(power, axis=0, keepdims=True)
    ) / np.maximum(np.std(power, axis=0, keepdims=True), 1e-12)
    eigenvalues = np.linalg.eigvalsh(
        np.corrcoef(standardized, rowvar=False)
    )
    common_mode_fraction = float(
        np.max(eigenvalues) / np.maximum(np.sum(eigenvalues), 1e-12)
    )
    return {
        "lag_seconds": lag_seconds,
        "cross_correlations": np.asarray(cross_correlations),
        "auto_correlations": np.asarray(auto_correlations),
        "delay_s": np.asarray(
            [value["delay_s"] for value in cross_metrics]
        ),
        "peak_correlation": np.asarray(
            [value["peak"] for value in cross_metrics]
        ),
        "peak_prominence": np.asarray(
            [value["prominence"] for value in cross_metrics]
        ),
        "peak_width_s": np.asarray(
            [value["width_s"] for value in cross_metrics]
        ),
        "peak_uniqueness": np.asarray(
            [value["uniqueness"] for value in cross_metrics]
        ),
        "common_mode_fraction": np.asarray(common_mode_fraction),
    }


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
    diagnostics = correlation_diagnostics(power, sample_interval_s)
    return (
        diagnostics["delay_s"],
        diagnostics["peak_correlation"],
    )


def fit_fading_window_diagnostics(
    power: np.ndarray,
    sample_interval_s: float,
) -> dict[str, np.ndarray | float | bool]:
    """
    Fit an evolving elliptical drifting-pattern correlation model.

    Fit all autocorrelation and cross-correlation curves, then report the
    fitted neutral wind and the objective checks used to accept or reject it.
    """
    power = np.asarray(power, dtype=np.float64)
    if power.ndim != 2 or power.shape[1] != len(CHANNEL_FIELDS):
        raise ValueError("power must have shape (time, three dipoles)")
    diagnostic = correlation_diagnostics(power, sample_interval_s)
    lag_seconds = diagnostic["lag_seconds"]
    if len(power) * sample_interval_s < 4 * MAX_LAG_S:
        raise ValueError("fading window is too short")

    cross_correlations = diagnostic["cross_correlations"]
    auto_correlations = diagnostic["auto_correlations"]
    observed = np.asarray(
        [np.mean(auto_correlations, axis=0), *cross_correlations]
    )
    zero_lag = int(np.argmin(np.abs(lag_seconds)))
    if 0 < zero_lag < observed.shape[1] - 1:
        # Receiver noise contributes only at zero lag in the autocorrelation.
        observed[0, zero_lag] = 0.5 * (
            observed[0, zero_lag - 1] + observed[0, zero_lag + 1]
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
        list(cross_correlations),
        lag_seconds,
    )

    def model_curves(parameters: np.ndarray) -> np.ndarray:
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
        return np.asarray(model_rows)

    residual_weight = np.sqrt(
        np.clip(np.abs(observed), 0.05, 1.0)
    )

    def model_and_residual(parameters: np.ndarray) -> np.ndarray:
        return ((model_curves(parameters) - observed) * residual_weight).ravel()

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
    model = model_curves(fit.x)
    model_peak_delays_s = np.asarray(
        [
            lag_seconds[int(np.nanargmax(value))]
            for value in model[1:]
        ],
    )
    model_peak_error_s = np.abs(
        model_peak_delays_s - diagnostic["delay_s"]
    )
    fit_rmse = float(
        np.sqrt(np.mean(model_and_residual(fit.x) ** 2))
    )
    velocity_jacobian = fit.jac[:, :2]
    jacobian_scale = np.maximum(
        np.linalg.norm(velocity_jacobian, axis=0, keepdims=True),
        1e-12,
    )
    singular_values = np.linalg.svd(
        velocity_jacobian / jacobian_scale,
        compute_uv=False,
    )
    fit_condition = float(
        singular_values[0] / max(singular_values[-1], 1e-12)
    )
    valid = (
        fit.success
        and np.all(
            diagnostic["peak_correlation"] >= MIN_BASELINE_CORRELATION
        )
        and np.all(
            diagnostic["peak_prominence"] >= MIN_PEAK_PROMINENCE
        )
        and np.all(
            diagnostic["peak_uniqueness"] >= MIN_PEAK_UNIQUENESS
        )
        and fit_rmse <= MAX_FIT_RMSE
        and pattern_speed <= MAX_PATTERN_SPEED_MS
        and np.all(np.abs(diagnostic["delay_s"]) < 0.95 * MAX_LAG_S)
        and np.all(model_peak_error_s <= MAX_MODEL_PEAK_ERROR_S)
        and float(diagnostic["common_mode_fraction"])
        <= MAX_COMMON_MODE_FRACTION
        and fit_condition <= MAX_SCALED_FIT_CONDITION
    )
    result = dict(diagnostic)
    result.update(
        {
            "valid": bool(valid),
            "zonal_wind_ms": (
                0.5 * float(pattern_velocity[0]) if valid else np.nan
            ),
            "meridional_wind_ms": (
                0.5 * float(pattern_velocity[1]) if valid else np.nan
            ),
            "candidate_zonal_wind_ms": 0.5 * float(pattern_velocity[0]),
            "candidate_meridional_wind_ms": 0.5 * float(pattern_velocity[1]),
            "fit_rmse": fit_rmse,
            "fit_condition": fit_condition,
            "model_curves": model,
            "model_peak_delays_s": model_peak_delays_s,
            "model_peak_error_s": model_peak_error_s,
        }
    )
    return result


def fit_fading_window(
    power: np.ndarray,
    sample_interval_s: float,
) -> tuple[float, float, float, float]:
    """Compatibility wrapper returning the accepted fit quantities."""
    result = fit_fading_window_diagnostics(power, sample_interval_s)
    return (
        float(result["zonal_wind_ms"]),
        float(result["meridional_wind_ms"]),
        float(np.median(result["peak_correlation"])),
        float(result["fit_rmse"]),
    )


def fit_robust_fading_window(
    power: np.ndarray,
    sample_interval_s: float,
) -> dict[str, np.ndarray | float | bool | int]:
    """Validate a 15-minute fit against independent five-minute subwindows."""
    full = fit_fading_window_diagnostics(power, sample_interval_s)
    subwindows = np.array_split(np.asarray(power), SUBWINDOW_COUNT, axis=0)
    subresults = [
        fit_fading_window_diagnostics(value, sample_interval_s)
        for value in subwindows
    ]
    subwind = np.asarray(
        [
            (
                value["zonal_wind_ms"],
                value["meridional_wind_ms"],
            )
            for value in subresults
        ],
        dtype=np.float64,
    )
    valid_subwindows = np.all(np.isfinite(subwind), axis=1)
    valid_values = subwind[valid_subwindows]
    valid_count = int(np.count_nonzero(valid_subwindows))
    component_spread = (
        np.ptp(valid_values, axis=0)
        if valid_count >= 2
        else np.full(2, np.nan)
    )
    bootstrap_sigma = np.full(2, np.nan)
    if valid_count >= 2:
        seed = int(
            np.round(
                np.sum(np.abs(valid_values)) * 1000.0
            )
        ) % (2**32)
        generator = np.random.default_rng(seed)
        bootstrap_index = generator.integers(
            0,
            valid_count,
            size=(BOOTSTRAP_REPLICATES, valid_count),
        )
        bootstrap_medians = np.median(
            valid_values[bootstrap_index],
            axis=1,
        )
        bootstrap_sigma = np.std(bootstrap_medians, axis=0, ddof=1)
    robust_valid = (
        bool(full["valid"])
        and valid_count >= MIN_VALID_SUBWINDOWS
        and np.all(
            component_spread <= MAX_SUBWINDOW_COMPONENT_SPREAD_MS
        )
        and np.all(
            bootstrap_sigma <= MAX_BOOTSTRAP_UNCERTAINTY_MS
        )
    )
    result = dict(full)
    result.update(
        {
            "valid": robust_valid,
            "zonal_wind_ms": (
                float(full["zonal_wind_ms"]) if robust_valid else np.nan
            ),
            "meridional_wind_ms": (
                float(full["meridional_wind_ms"]) if robust_valid else np.nan
            ),
            "valid_subwindow_count": valid_count,
            "subwindow_component_spread_ms": component_spread,
            "bootstrap_uncertainty_ms": bootstrap_sigma,
        }
    )
    return result


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
