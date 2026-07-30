#!/usr/bin/env python3
"""Generate dense fitted-Doppler RTIs for the latest 30 minutes."""

from __future__ import annotations

import datetime as dt
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
from plot_dense_doppler import dense_doppler, fit_sinusoid_bank


WINDOW_S = 30 * 60
DISPLAY_LIMIT_MS = 200.0
DISPLAY_RANGE_MAX_KM = 300.0
PLOT_DIR = Path("/data2/plots/monitor")
COMBINED_PLOT = PLOT_DIR / "latest_snr_doppler_30m_0_300.png"
CHANNEL_HEALTH_PLOT = PLOT_DIR / "latest_channel_snr_30m_0_300.png"
POSITIONS_PLOT = PLOT_DIR / "latest_dense_positions_5m.png"
DATA_FILE = PLOT_DIR / "latest_doppler_30m.npz"
SNR_STATE_FILE = PLOT_DIR / "rti_30m_2s_state.npz"
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


def load_channel_power(
    reader: drf.DigitalMetadataReader,
    start_unix: float,
    end_unix: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return two-second power for every receiver channel and range gate."""
    records = []
    start_us = int(start_unix * 1e6)
    end_us = int(end_unix * 1e6)
    for chunk_start in np.arange(start_us, end_us, int(60e6)):
        chunk_end = min(chunk_start + int(60e6), end_us)
        metadata = reader.read(int(chunk_start), int(chunk_end) - 1)
        for key in sorted(metadata):
            record = metadata[key]
            fields = [field for field, _ in CHANNEL_HEALTH_FIELDS]
            if not all(field in record for field in (*fields, "rvec")):
                continue
            ranges = np.asarray(record["rvec"], dtype=np.float64)
            mask = (ranges >= 0.0) & (ranges <= DISPLAY_RANGE_MAX_KM)
            if not np.any(mask):
                continue
            power = np.asarray(
                [
                    np.mean(
                        np.abs(np.asarray(record[field])[:, mask]) ** 2,
                        axis=0,
                    )
                    for field in fields
                ],
                dtype=np.float32,
            )
            records.append(
                (float(key) / 1e6 + 1.0, ranges[mask], power)
            )

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
        "Ramfjordmoen MF radar · receiver-channel health · latest 30 minutes",
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
    doppler_times: np.ndarray,
    doppler_ranges: np.ndarray,
    velocity_ms: np.ndarray,
    start: dt.datetime,
    end: dt.datetime,
) -> None:
    """Plot synchronized two-second SNR and fitted-Doppler RTI panels."""
    with np.load(SNR_STATE_FILE) as state:
        snr_times = state["times"].astype(np.int64, copy=False)
        snr_power = state["power"].astype(np.float32, copy=False)

    start_unix = start.timestamp()
    end_unix = end.timestamp()
    time_mask = (snr_times >= start_unix) & (snr_times <= end_unix)
    snr_times = snr_times[time_mask]
    snr_power = snr_power[time_mask]
    if len(snr_times) < 2:
        raise RuntimeError("two-second SNR state does not overlap Doppler window")

    snr_ranges = monitor.range_vector_km(snr_power.shape[1])
    noise_mask = (
        (snr_ranges >= monitor.THIRTY_MINUTE_NOISE_RANGE_KM[0])
        & (snr_ranges <= monitor.THIRTY_MINUTE_NOISE_RANGE_KM[1])
    )
    display_mask = (
        (snr_ranges >= 0.0)
        & (snr_ranges <= DISPLAY_RANGE_MAX_KM)
    )
    background_power = max(
        float(np.nanmean(snr_power[:, noise_mask])),
        1e-20,
    )
    power_ratio_db = 10.0 * np.log10(
        np.maximum(
            snr_power[:, display_mask] / background_power,
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
        power_ratio_db.T,
        shading="nearest",
        cmap="plasma",
        vmin=monitor.THIRTY_MINUTE_POWER_RATIO_MIN_DB,
        vmax=20.0,
    )
    axes[0].set_title(
        "Power / 30–50 km interval background · two-second cadence"
    )
    snr_colorbar = figure.colorbar(snr_image, ax=axes[0], pad=0.01)
    snr_colorbar.set_label("Power / background (dB)", color="#b7c5d9")
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
        "Unfiltered one-second complex-sinusoid fit · every cell"
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
        "Ramfjordmoen MF radar · latest 30 minutes · "
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

    times, ranges, velocity, fit_snr = dense_doppler(
        reader,
        start_unix,
        end_unix,
        (0.0, DISPLAY_RANGE_MAX_KM),
    )
    channel_times, channel_ranges, channel_power = load_channel_power(
        reader,
        start_unix,
        end_unix,
    )
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    plot_channel_health(
        channel_times,
        channel_ranges,
        channel_power,
        start,
        end,
    )
    plot_combined_rti(
        times,
        ranges,
        velocity,
        start,
        end,
    )
    np.savez_compressed(
        DATA_FILE,
        time_unix=times,
        range_km=ranges,
        velocity_ms=velocity,
        sinusoid_snr=fit_snr,
    )
    print(f"Four-channel SNR health: {CHANNEL_HEALTH_PLOT}")
    print(f"Combined SNR/Doppler RTI: {COMBINED_PLOT}")
    print(f"Dense fit cells: {velocity.size}")


if __name__ == "__main__":
    main()
