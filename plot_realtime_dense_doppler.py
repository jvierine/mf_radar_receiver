#!/usr/bin/env python3
"""Generate dense fitted-Doppler RTIs for the latest 15 minutes."""

from __future__ import annotations

import datetime as dt
import h5py
import os
from pathlib import Path

import digital_rf as drf
import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

import image_help as ih
import mf_conf as mc
import plot_monitor_rti as monitor
import wind_estimates_3ch_2days as winds
from plot_dense_doppler import fit_common_sinusoid_fft, fit_sinusoid_bank


WINDOW_S = 15 * 60
FIT_DURATION_S = 10
SOURCE_RECORD_S = 2
DISPLAY_LIMIT_MS = 100.0
DISPLAY_RANGE_MAX_KM = 300.0
NARROWBAND_RATIO_MIN_DB = -20.0
PLOT_DIR = Path("/data2/plots/monitor")
COMBINED_PLOT = PLOT_DIR / "latest_snr_doppler_15m_0_300.png"
CHANNEL_HEALTH_PLOT = PLOT_DIR / "latest_channel_snr_15m_0_300.png"
POSITIONS_PLOT = PLOT_DIR / "latest_dense_positions_5m.png"
DATA_FILE = PLOT_DIR / "latest_doppler_15m.h5"
POSITION_WINDOW_S = 5 * 60
POSITION_FIT_SNR_MIN = 5.0
POSITION_COHERENCE_MIN = 0.8
POSITION_PHASE_CLOSURE_MAX = 0.35
POSITION_PHASE_RESIDUAL_MAX = 0.35
CHANNEL_HEALTH_FIELDS = (
    ("rti1", "Channel 1 · MF3 dipole"),
    ("rti2", "Channel 2 · loop"),
    ("rti3", "Channel 3 · MF1 dipole"),
    ("rti4", "Channel 4 · MF2 dipole"),
)
DOPPLER_FIELDS = ("rti1", "rti3", "rti4")
DOPPLER_CHANNEL_INDICES = np.asarray((0, 2, 3))


def iter_fit_groups(
    reader: drf.DigitalMetadataReader,
    start_unix: float,
    end_unix: float,
):
    """Yield complete, non-overlapping groups spanning ten seconds."""
    first_start = int(np.ceil(start_unix / FIT_DURATION_S) * FIT_DURATION_S)
    final_start = int(np.floor(end_unix / FIT_DURATION_S) * FIT_DURATION_S)
    for group_start in range(first_start, final_start, FIT_DURATION_S):
        expected_keys = [
            int((group_start + offset) * 1e6)
            for offset in range(0, FIT_DURATION_S, SOURCE_RECORD_S)
        ]
        metadata = reader.read(
            expected_keys[0],
            int((group_start + FIT_DURATION_S) * 1e6) - 1,
        )
        if not all(key in metadata for key in expected_keys):
            continue
        records = [metadata[key] for key in expected_keys]
        if not all(
            field in record
            for record in records
            for field in ("rvec", "tvec")
        ):
            continue
        reference_range = np.asarray(records[0]["rvec"], dtype=np.float64)
        if not all(
            np.array_equal(
                reference_range,
                np.asarray(record["rvec"], dtype=np.float64),
            )
            for record in records[1:]
        ):
            continue
        yield group_start + FIT_DURATION_S / 2.0, expected_keys, records


def load_channel_power(
    reader: drf.DigitalMetadataReader,
    start_unix: float,
    end_unix: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return full-ten-second power for every channel and range gate."""
    records = []
    fields = [field for field, _ in CHANNEL_HEALTH_FIELDS]
    for center_time, _, group in iter_fit_groups(
        reader,
        start_unix,
        end_unix,
    ):
        if not all(field in record for record in group for field in fields):
            continue
        ranges = np.asarray(group[0]["rvec"], dtype=np.float64)
        mask = (ranges >= 0.0) & (ranges <= DISPLAY_RANGE_MAX_KM)
        if not np.any(mask):
            continue
        power = np.asarray(
            [
                np.mean(
                    np.abs(
                        np.concatenate(
                            [
                                np.asarray(record[field])[:, mask]
                                for record in group
                            ],
                            axis=0,
                        )
                    )
                    ** 2,
                    axis=0,
                )
                for field in fields
            ],
            dtype=np.float32,
        )
        records.append((center_time, ranges[mask], power))

    if not records:
        raise RuntimeError("no four-channel metadata in requested interval")
    range_km = max((row[1] for row in records), key=len)
    power = np.full(
        (len(records), len(CHANNEL_HEALTH_FIELDS), len(range_km)),
        np.nan,
        dtype=np.float32,
    )
    for row_index, (_, record_range, record_power) in enumerate(records):
        indices = np.searchsorted(range_km, record_range)
        if (
            np.any(indices >= len(range_km))
            or not np.allclose(range_km[indices], record_range)
        ):
            raise ValueError("incompatible range vectors within interval")
        power[row_index][:, indices] = record_power
    return (
        np.asarray([row[0] for row in records], dtype=np.float64),
        range_km,
        power,
    )


def fit_ten_second_doppler(
    reader: drf.DigitalMetadataReader,
    start_unix: float,
    end_unix: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fit one sinusoid to every channel/range cell over each full 10 s."""
    rows = []
    fields = DOPPLER_FIELDS
    for center_time, keys, group in iter_fit_groups(
        reader,
        start_unix,
        end_unix,
    ):
        if not all(field in record for record in group for field in fields):
            continue
        ranges = np.asarray(group[0]["rvec"], dtype=np.float64)
        mask = (ranges >= 0.0) & (ranges <= DISPLAY_RANGE_MAX_KM)
        if not np.any(mask):
            continue
        selected_range = ranges[mask]
        fit_times = np.concatenate(
            [
                float(key) / 1e6
                + np.asarray(record["tvec"], dtype=np.float64)
                for key, record in zip(keys, group)
            ]
        )
        fit_times -= fit_times[0]
        channel_voltage = np.asarray(
            [
                np.concatenate(
                    [
                        np.asarray(record[field])[:, mask]
                        for record in group
                    ],
                    axis=0,
                )
                for field in fields
            ],
            dtype=np.complex128,
        )
        frequency, snr, fitted_amplitude = fit_common_sinusoid_fft(
            fit_times,
            channel_voltage.transpose(1, 0, 2),
            winds.MAX_FIT_DOPPLER_HZ,
        )
        rows.append(
            (
                center_time,
                selected_range,
                winds.DOPPLER_SIGN
                * winds.DOPPLER_TO_MS
                * frequency,
                snr,
                np.abs(fitted_amplitude) ** 2,
            )
        )

    if not rows:
        raise RuntimeError("no complete ten-second metadata groups")
    range_km = max((row[1] for row in rows), key=len)
    velocity = np.full((len(rows), len(range_km)), np.nan)
    fit_snr = np.full_like(velocity, np.nan)
    fitted_power = np.full(
        (len(rows), len(DOPPLER_FIELDS), len(range_km)),
        np.nan,
    )
    for row_index, (
        _,
        record_range,
        record_velocity,
        record_snr,
        record_fitted_power,
    ) in enumerate(rows):
        indices = np.searchsorted(range_km, record_range)
        if (
            np.any(indices >= len(range_km))
            or not np.allclose(range_km[indices], record_range)
        ):
            raise ValueError("incompatible range vectors within interval")
        velocity[row_index, indices] = record_velocity
        fit_snr[row_index, indices] = record_snr
        fitted_power[row_index][:, indices] = record_fitted_power
    return (
        np.asarray([row[0] for row in rows], dtype=np.float64),
        range_km,
        velocity,
        fit_snr,
        fitted_power,
    )


def plot_channel_health(
    times_unix: np.ndarray,
    range_km: np.ndarray,
    power: np.ndarray,
    start: dt.datetime,
    end: dt.datetime,
) -> None:
    """Plot separately noise-referenced power for all four receiver channels."""
    noise_mask = (
        (range_km >= monitor.THIRTY_MINUTE_NOISE_RANGE_KM[0])
        & (range_km <= monitor.THIRTY_MINUTE_NOISE_RANGE_KM[1])
    )
    background_power = np.nanmean(
        power[:, :, noise_mask],
        axis=(0, 2),
    )
    power_ratio_db = 10.0 * np.log10(
        np.maximum(
            power / np.maximum(background_power[None, :, None], 1e-20),
            1e-20,
        )
    )
    time_values = [
        dt.datetime.fromtimestamp(value, tz=dt.timezone.utc)
        for value in times_unix
    ]
    figure, axes = plt.subplots(
        4,
        1,
        figsize=(15, 12),
        sharex=True,
        sharey=True,
        facecolor="#070b14",
        constrained_layout=True,
    )
    image = None
    for channel_index, (axis, (_, label)) in enumerate(
        zip(axes, CHANNEL_HEALTH_FIELDS)
    ):
        axis.set_facecolor("#050810")
        image = axis.pcolormesh(
            time_values,
            range_km,
            power_ratio_db[:, channel_index, :].T,
            shading="nearest",
            cmap="plasma",
            vmin=monitor.THIRTY_MINUTE_POWER_RATIO_MIN_DB,
            vmax=20.0,
        )
        noise_db = 10.0 * np.log10(
            max(float(background_power[channel_index]), 1e-20)
        )
        axis.set_title(
            f"{label} · background {noise_db:.1f} dB (arb. power)",
            color="#edf4ff",
        )
        axis.set_xlim(start, end)
        axis.set_ylim(0.0, DISPLAY_RANGE_MAX_KM)
        axis.set_ylabel("Round-trip range (km)", color="#b7c5d9")
        axis.tick_params(colors="#8fa1ba")
        for spine in axis.spines.values():
            spine.set_color("#314563")
    axes[-1].set_xlabel("Time (UTC)", color="#b7c5d9")
    axes[-1].xaxis.set_major_locator(mdates.MinuteLocator(interval=5))
    axes[-1].xaxis.set_major_formatter(
        mdates.DateFormatter("%H:%M", tz=dt.timezone.utc)
    )
    colorbar = figure.colorbar(image, ax=axes, pad=0.01)
    colorbar.set_label(
        "Power / channel background (dB)",
        color="#b7c5d9",
    )
    colorbar.ax.tick_params(colors="#8fa1ba")
    figure.suptitle(
        "Ramfjordmoen MF radar · receiver-channel health · latest 15 minutes",
        color="#edf4ff",
        fontsize=18,
        weight="semibold",
    )
    temporary = CHANNEL_HEALTH_PLOT.with_suffix(".tmp.png")
    figure.savefig(
        temporary,
        dpi=150,
        facecolor=figure.get_facecolor(),
    )
    plt.close(figure)
    os.replace(temporary, CHANNEL_HEALTH_PLOT)


def plot_combined_rti(
    snr_times: np.ndarray,
    snr_ranges: np.ndarray,
    channel_power: np.ndarray,
    fitted_channel_power: np.ndarray,
    doppler_times: np.ndarray,
    doppler_ranges: np.ndarray,
    velocity_ms: np.ndarray,
    start: dt.datetime,
    end: dt.datetime,
) -> None:
    """Plot synchronized full-ten-second SNR and fitted-Doppler panels."""
    if len(snr_times) < 2:
        raise RuntimeError("fewer than two ten-second SNR estimates")

    noise_mask = (
        (snr_ranges >= monitor.THIRTY_MINUTE_NOISE_RANGE_KM[0])
        & (snr_ranges <= monitor.THIRTY_MINUTE_NOISE_RANGE_KM[1])
    )
    display_mask = (
        (snr_ranges >= 0.0)
        & (snr_ranges <= DISPLAY_RANGE_MAX_KM)
    )
    broadband_noise_power = np.maximum(
        np.nanmean(channel_power[:, :, noise_mask], axis=(0, 2)),
        1e-20,
    )
    doppler_noise_power = broadband_noise_power[DOPPLER_CHANNEL_INDICES]
    narrowband_ratio_db = 10.0 * np.log10(
        np.maximum(
            np.sum(
                fitted_channel_power[:, :, display_mask]
                / doppler_noise_power[None, :, None],
                axis=1,
            ),
            1e-20,
        )
    )

    figure, axes = plt.subplots(
        2,
        1,
        figsize=(15, 9),
        sharex=True,
        sharey=True,
        facecolor="#070b14",
        constrained_layout=True,
    )
    for axis in axes:
        axis.set_facecolor("#050810")
        axis.tick_params(colors="#8fa1ba")
        axis.yaxis.label.set_color("#b7c5d9")
        axis.title.set_color("#edf4ff")
        for spine in axis.spines.values():
            spine.set_color("#314563")

    snr_image = axes[0].pcolormesh(
        [
            dt.datetime.fromtimestamp(value, tz=dt.timezone.utc)
            for value in snr_times
        ],
        snr_ranges[display_mask],
        narrowband_ratio_db.T,
        shading="nearest",
        cmap="plasma",
        vmin=NARROWBAND_RATIO_MIN_DB,
        vmax=20.0,
    )
    axes[0].set_title(
        "Joint fitted narrowband power / processed broadband noise"
    )
    snr_colorbar = figure.colorbar(snr_image, ax=axes[0], pad=0.01)
    snr_colorbar.set_label(
        "Σ fitted power / broadband noise (dB)",
        color="#b7c5d9",
    )
    snr_colorbar.ax.tick_params(colors="#8fa1ba")

    doppler_image = axes[1].pcolormesh(
        [
            dt.datetime.fromtimestamp(value, tz=dt.timezone.utc)
            for value in doppler_times
        ],
        doppler_ranges,
        velocity_ms.T,
        shading="nearest",
        cmap="seismic",
        vmin=-DISPLAY_LIMIT_MS,
        vmax=DISPLAY_LIMIT_MS,
    )
    axes[1].set_title(
        "Joint three-dipole 10-second complex-sinusoid fit · every cell"
    )
    doppler_colorbar = figure.colorbar(
        doppler_image,
        ax=axes[1],
        pad=0.01,
    )
    doppler_colorbar.set_label(
        "Monostatic Doppler velocity (m/s)",
        color="#b7c5d9",
    )
    doppler_colorbar.ax.tick_params(colors="#8fa1ba")

    for axis in axes:
        axis.set_xlim(start, end)
        axis.set_ylim(0.0, DISPLAY_RANGE_MAX_KM)
        axis.set_ylabel("Round-trip range (km)")
    axes[1].set_xlabel("Time (UTC)", color="#b7c5d9")
    axes[1].xaxis.set_major_locator(mdates.MinuteLocator(interval=5))
    axes[1].xaxis.set_major_formatter(
        mdates.DateFormatter("%H:%M", tz=dt.timezone.utc)
    )
    figure.suptitle(
        "Ramfjordmoen MF radar · latest 15 minutes · "
        "SNR and fitted Doppler",
        color="#edf4ff",
        fontsize=18,
        weight="semibold",
    )
    temporary = COMBINED_PLOT.with_suffix(".tmp.png")
    figure.savefig(
        temporary,
        dpi=150,
        facecolor=figure.get_facecolor(),
    )
    plt.close(figure)
    os.replace(temporary, COMBINED_PLOT)


def plot_recent_positions(
    reader: drf.DigitalMetadataReader,
    end_unix: float,
) -> int:
    """Image the latest five minutes directly from dense per-cell fits."""
    phasecal = winds.load_phasecal()
    coordinates = winds.antenna_coordinates_ecef()
    position_differences = [
        coordinates[0, :] - coordinates[2, :],
        coordinates[0, :] - coordinates[3, :],
        coordinates[2, :] - coordinates[3, :],
    ]
    position_start = end_unix - POSITION_WINDOW_S
    metadata = reader.read(
        int(position_start * 1e6),
        int(end_unix * 1e6) - 1,
    )
    position_rows = []
    for key in sorted(metadata):
        record = metadata[key]
        if not all(
            field in record
            for field in ("rti1", "rti3", "rti4", "rvec", "tvec")
        ):
            continue
        ranges = np.asarray(record["rvec"], dtype=np.float64)
        range_mask = (ranges >= 70.0) & (ranges <= 200.0)
        ranges = ranges[range_mask]
        if not len(ranges):
            continue
        times = np.asarray(record["tvec"], dtype=np.float64)
        channel_voltage = np.asarray(
            [
                np.asarray(record[field])[:, range_mask]
                * np.exp(1j * phasecal[channel_index])
                for field, channel_index in (
                    ("rti1", 0),
                    ("rti3", 2),
                    ("rti4", 3),
                )
            ],
            dtype=np.complex128,
        )
        sample_interval = float(np.median(np.diff(times)))
        fit_samples = min(
            len(times),
            max(
                2,
                int(
                    round(
                        winds.DOPPLER_FIT_DURATION_S / sample_interval
                    )
                ),
            ),
        )
        fit_start = (len(times) - fit_samples) // 2
        fit_stop = fit_start + fit_samples
        fit_times = times[fit_start:fit_stop] - times[fit_start]
        fit_voltage = channel_voltage[:, fit_start:fit_stop, :]
        scale = np.sqrt(np.mean(np.abs(fit_voltage) ** 2, axis=1))
        valid_scale = np.isfinite(scale) & (scale > 0)
        safe_scale = np.where(valid_scale, scale, 1.0)
        normalized = fit_voltage / safe_scale[:, None, :]

        gate_count = len(ranges)
        bank = normalized.transpose(1, 0, 2).reshape(
            fit_samples,
            3 * gate_count,
        )
        frequency, fit_snr = fit_sinusoid_bank(fit_times, bank)
        frequency = frequency.reshape(3, gate_count)
        fit_snr = fit_snr.reshape(3, gate_count)
        best_channel = np.argmax(
            np.where(np.isfinite(fit_snr), fit_snr, -np.inf),
            axis=0,
        )
        gate = np.arange(gate_count)
        best_frequency = frequency[best_channel, gate]
        best_snr = fit_snr[best_channel, gate]

        xc13 = np.mean(normalized[0] * np.conj(normalized[1]), axis=0)
        xc14 = np.mean(normalized[0] * np.conj(normalized[2]), axis=0)
        xc34 = np.mean(normalized[1] * np.conj(normalized[2]), axis=0)
        cross = np.asarray([xc13, xc14, xc34])
        coherence = np.mean(np.abs(cross), axis=0)
        closure = np.abs(np.angle(xc13 * np.conj(xc14) * xc34))
        candidate_gates = np.flatnonzero(
            np.all(valid_scale, axis=0)
            & np.isfinite(best_frequency)
            & (best_snr >= POSITION_FIT_SNR_MIN)
            & (coherence >= POSITION_COHERENCE_MIN)
            & (closure <= POSITION_PHASE_CLOSURE_MAX)
        )
        for range_index in candidate_gates:
            candidates = ih.aoa_candidates(
                cross[:, range_index],
                position_differences,
                wavelength=mc.wavelength,
            )
            accepted = []
            for east, north, up, phase_residual, match in candidates:
                altitude = ranges[range_index] * up
                if not 70.0 <= altitude <= 150.0:
                    continue
                if phase_residual > POSITION_PHASE_RESIDUAL_MAX:
                    continue
                cost = (
                    (phase_residual / POSITION_PHASE_RESIDUAL_MAX) ** 2
                    + (1.0 - match)
                )
                accepted.append((cost, east, north, up))
            if not accepted:
                continue
            _, east, north, up = min(accepted, key=lambda row: row[0])
            radial_velocity = (
                winds.DOPPLER_SIGN
                * winds.DOPPLER_TO_MS
                * best_frequency[range_index]
            )
            position_rows.append(
                (
                    float(key) / 1e6 + 1.0,
                    ranges[range_index] * east,
                    ranges[range_index] * north,
                    ranges[range_index] * up,
                    radial_velocity,
                    best_snr[range_index],
                    coherence[range_index],
                    ranges[range_index],
                )
            )
    detections = np.asarray(position_rows, dtype=np.float64)
    if detections.size == 0:
        detections = np.empty((0, 8), dtype=np.float64)

    figure, axis = plt.subplots(
        figsize=(8, 8),
        facecolor="#070b14",
        constrained_layout=True,
    )
    axis.set_facecolor("#050810")
    axis.set_xlim(-200, 200)
    axis.set_ylim(-200, 200)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("East–west position (km)", color="#b7c5d9")
    axis.set_ylabel("North–south position (km)", color="#b7c5d9")
    axis.tick_params(colors="#8fa1ba")
    axis.grid(color="#ffffff", alpha=0.07)
    for spine in axis.spines.values():
        spine.set_color("#314563")
    axis.scatter(
        [0],
        [0],
        marker="+",
        s=110,
        linewidths=1.5,
        color="#edf4ff",
    )

    if len(detections):
        velocity_limit = float(
            np.clip(
                np.ceil(np.nanpercentile(np.abs(detections[:, 4]), 98)),
                20,
                300,
            )
        )
        positions = axis.scatter(
            detections[:, 1],
            detections[:, 2],
            c=detections[:, 4],
            s=18,
            cmap="seismic",
            vmin=-velocity_limit,
            vmax=velocity_limit,
            linewidths=0,
            alpha=0.9,
            rasterized=True,
        )
        colorbar = figure.colorbar(positions, ax=axis, pad=0.02)
        colorbar.set_label(
            "Monostatic Doppler velocity (m/s)",
            color="#b7c5d9",
        )
        colorbar.ax.tick_params(colors="#8fa1ba")
    else:
        axis.text(
            0.5,
            0.5,
            "No position-quality echoes in the latest five minutes",
            transform=axis.transAxes,
            ha="center",
            va="center",
            color="#b7c5d9",
        )
    axis.set_title(
        f"Interferometric positions · 5 min · N={len(detections)}",
        color="#edf4ff",
        fontsize=14,
        weight="semibold",
    )
    temporary = POSITIONS_PLOT.with_suffix(".tmp.png")
    figure.savefig(
        temporary,
        dpi=170,
        facecolor=figure.get_facecolor(),
    )
    plt.close(figure)
    os.replace(temporary, POSITIONS_PLOT)
    return len(detections)


def main() -> None:
    reader = drf.DigitalMetadataReader(mc.xc_dir)
    bounds = reader.get_bounds()
    end_unix = float(bounds[1]) / 1e6 + 2.0
    start_unix = end_unix - WINDOW_S
    start = dt.datetime.fromtimestamp(start_unix, tz=dt.timezone.utc)
    end = dt.datetime.fromtimestamp(end_unix, tz=dt.timezone.utc)

    times, ranges, velocity, fit_snr, fitted_channel_power = (
        fit_ten_second_doppler(
        reader,
        start_unix,
        end_unix,
        )
    )
    channel_times, channel_ranges, channel_power = load_channel_power(
        reader,
        start_unix,
        end_unix,
    )
    if not np.array_equal(times, channel_times):
        raise RuntimeError("SNR and Doppler ten-second timestamps differ")
    if not np.array_equal(ranges, channel_ranges):
        raise RuntimeError("SNR and Doppler range vectors differ")
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    plot_channel_health(
        channel_times,
        channel_ranges,
        channel_power,
        start,
        end,
    )
    plot_combined_rti(
        channel_times,
        channel_ranges,
        channel_power,
        fitted_channel_power,
        times,
        ranges,
        velocity,
        start,
        end,
    )
    temporary_data = DATA_FILE.with_suffix(DATA_FILE.suffix + ".tmp")
    with h5py.File(temporary_data, "w") as handle:
        handle.attrs["fit_duration_seconds"] = FIT_DURATION_S
        handle.attrs["output_cadence_seconds"] = FIT_DURATION_S
        handle.attrs["doppler_combination"] = (
            "joint_common_frequency_independent_complex_amplitudes"
        )
        handle.attrs["doppler_channels"] = "ch1,ch3,ch4"
        handle.attrs["excluded_doppler_channels"] = (
            "ch2 loop: lower SNR and stronger RFI"
        )
        handle.attrs["narrowband_power"] = (
            "squared_complex_amplitude_of_joint_common_frequency_fit"
        )
        handle.attrs["broadband_noise_power"] = (
            "mean_processed_voltage_power_over_30_to_50_km_and_15_minutes"
        )
        handle.create_dataset("time_unix", data=times)
        handle.create_dataset("range_km", data=ranges)
        handle.create_dataset(
            "velocity_ms",
            data=velocity,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
        handle.create_dataset(
            "sinusoid_snr",
            data=fit_snr,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
        handle.create_dataset(
            "fitted_narrowband_power_by_channel",
            data=fitted_channel_power,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
        noise_mask = (
            (channel_ranges >= monitor.THIRTY_MINUTE_NOISE_RANGE_KM[0])
            & (channel_ranges <= monitor.THIRTY_MINUTE_NOISE_RANGE_KM[1])
        )
        broadband_noise_power = np.nanmean(
            channel_power[:, :, noise_mask],
            axis=(0, 2),
        )
        doppler_noise_power = broadband_noise_power[DOPPLER_CHANNEL_INDICES]
        handle.create_dataset(
            "processed_broadband_noise_power_by_channel",
            data=doppler_noise_power,
        )
        handle.create_dataset(
            "narrowband_to_broadband_noise_ratio",
            data=np.sum(
                fitted_channel_power
                / np.maximum(
                    doppler_noise_power[None, :, None],
                    1e-20,
                ),
                axis=1,
            ),
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
    os.replace(temporary_data, DATA_FILE)
    print(f"Four-channel SNR health: {CHANNEL_HEALTH_PLOT}")
    print(f"Combined SNR/Doppler RTI: {COMBINED_PLOT}")
    print(f"Dense fit cells: {velocity.size}")


if __name__ == "__main__":
    main()
