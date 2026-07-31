#!/usr/bin/env python3
"""Generate dense fitted-Doppler RTIs for the latest 15 minutes."""

from __future__ import annotations

import datetime as dt
import h5py
import os
from pathlib import Path
import subprocess
import sys

import digital_rf as drf
import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

import fading_wind
import mf_conf as mc
import plot_fading_wind_history as fading_history
import plot_monitor_rti as monitor
import radar_common as common
from plot_dense_doppler import fit_common_sinusoid_fft


WINDOW_S = 15 * 60
FIT_DURATION_S = 10
SOURCE_RECORD_S = 2
DISPLAY_LIMIT_MS = 150.0
DISPLAY_RANGE_MAX_KM = 300.0
NARROWBAND_RATIO_MIN_DB = -20.0
PLOT_DIR = Path("/data2/plots/monitor")
COMBINED_PLOT = PLOT_DIR / "latest_snr_doppler_15m_0_300.png"
CHANNEL_HEALTH_PLOT = PLOT_DIR / "latest_channel_snr_15m_0_300.png"
DATA_FILE = PLOT_DIR / "latest_doppler_15m.h5"
FADING_DIAGNOSTIC_PLOT = (
    PLOT_DIR / "latest_fading_correlation_diagnostics_15m.png"
)
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
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Fit one common Doppler to all three dipoles in every range cell."""
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
        frequency, fit_snr, fitted_amplitude = fit_common_sinusoid_fft(
            fit_times,
            channel_voltage.transpose(1, 0, 2),
            common.MAX_FIT_DOPPLER_HZ,
        )
        rows.append(
            (
                center_time,
                selected_range,
                frequency,
                common.DOPPLER_SIGN * common.DOPPLER_TO_MS * frequency,
                fit_snr,
                np.abs(fitted_amplitude) ** 2,
            )
        )

    if not rows:
        raise RuntimeError("no complete ten-second metadata groups")
    range_km = max((row[1] for row in rows), key=len)
    frequency = np.full((len(rows), len(range_km)), np.nan)
    velocity = np.full((len(rows), len(range_km)), np.nan)
    fit_snr = np.full_like(velocity, np.nan)
    fitted_power = np.full(
        (len(rows), len(DOPPLER_FIELDS), len(range_km)),
        np.nan,
    )
    for row_index, (
        _,
        record_range,
        record_frequency,
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
        frequency[row_index, indices] = record_frequency
        velocity[row_index, indices] = record_velocity
        fit_snr[row_index, indices] = record_snr
        fitted_power[row_index][:, indices] = record_fitted_power
    return (
        np.asarray([row[0] for row in rows], dtype=np.float64),
        range_km,
        frequency,
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
    power_times: np.ndarray,
    power_ranges: np.ndarray,
    narrowband_power_ratio: np.ndarray,
    doppler_times: np.ndarray,
    doppler_ranges: np.ndarray,
    velocity_ms: np.ndarray,
    wind_times: np.ndarray,
    zonal_wind_ms: np.ndarray,
    meridional_wind_ms: np.ndarray,
    start: dt.datetime,
    end: dt.datetime,
) -> None:
    """Plot ungated RTIs and quality-controlled fading-wind estimates."""
    if len(power_times) < 2:
        raise RuntimeError("fewer than two ten-second SNR estimates")

    display_mask = (
        (power_ranges >= 0.0)
        & (power_ranges <= DISPLAY_RANGE_MAX_KM)
    )
    narrowband_ratio_db = 10.0 * np.log10(
        np.maximum(
            narrowband_power_ratio[:, display_mask],
            1e-20,
        )
    )

    figure, axes = plt.subplots(
        4,
        1,
        figsize=(15, 16),
        sharex=True,
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
            for value in power_times
        ],
        power_ranges[display_mask],
        narrowband_ratio_db.T,
        shading="nearest",
        cmap="plasma",
        vmin=NARROWBAND_RATIO_MIN_DB,
        vmax=20.0,
    )
    axes[0].set_title(
        "Three-dipole fitted narrowband power / processed broadband noise"
    )
    snr_colorbar = figure.colorbar(snr_image, ax=axes[0], pad=0.01)
    snr_colorbar.set_label(
        "Σ fitted dipole power / broadband noise (dB)",
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
        "Three-dipole common-frequency 10-second fit · every cell retained"
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

    wind_datetimes = [
        dt.datetime.fromtimestamp(value, tz=dt.timezone.utc)
        for value in wind_times
    ]
    zonal_image = axes[2].pcolormesh(
        wind_datetimes,
        fading_wind.ALTITUDE_KM,
        zonal_wind_ms.T,
        shading="nearest",
        cmap="seismic",
        vmin=-DISPLAY_LIMIT_MS,
        vmax=DISPLAY_LIMIT_MS,
    )
    axes[2].set_title(
        "Zonal neutral wind · east positive · five-minute fading correlation"
    )
    zonal_colorbar = figure.colorbar(zonal_image, ax=axes[2], pad=0.01)
    zonal_colorbar.set_label("Zonal wind (m/s)", color="#b7c5d9")
    zonal_colorbar.ax.tick_params(colors="#8fa1ba")

    meridional_image = axes[3].pcolormesh(
        wind_datetimes,
        fading_wind.ALTITUDE_KM,
        meridional_wind_ms.T,
        shading="nearest",
        cmap="seismic",
        vmin=-DISPLAY_LIMIT_MS,
        vmax=DISPLAY_LIMIT_MS,
    )
    axes[3].set_title(
        "Meridional neutral wind · north positive · five-minute fading correlation"
    )
    meridional_colorbar = figure.colorbar(
        meridional_image,
        ax=axes[3],
        pad=0.01,
    )
    meridional_colorbar.set_label(
        "Meridional wind (m/s)",
        color="#b7c5d9",
    )
    meridional_colorbar.ax.tick_params(colors="#8fa1ba")

    for axis in axes[:2]:
        axis.set_xlim(start, end)
        axis.set_ylim(0.0, DISPLAY_RANGE_MAX_KM)
        axis.set_ylabel("Round-trip range (km)")
    for axis in axes[2:]:
        axis.set_xlim(start, end)
        axis.set_ylim(
            fading_wind.ALTITUDE_KM[0]
            - fading_wind.ALTITUDE_HALF_WIDTH_KM,
            fading_wind.ALTITUDE_KM[-1]
            + fading_wind.ALTITUDE_HALF_WIDTH_KM,
        )
        axis.set_ylabel("Altitude (km)")
    axes[-1].set_xlabel("Time (UTC)", color="#b7c5d9")
    axes[-1].xaxis.set_major_locator(mdates.MinuteLocator(interval=5))
    axes[-1].xaxis.set_major_formatter(
        mdates.DateFormatter("%H:%M", tz=dt.timezone.utc)
    )
    figure.suptitle(
        "Ramfjordmoen MF radar · latest 15 minutes · "
        "power, Doppler, and fading winds",
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


def plot_fading_correlation_diagnostics(
    wind_times: np.ndarray,
    baseline_delay_s: np.ndarray,
    baseline_peak_correlation: np.ndarray,
    start: dt.datetime,
    end: dt.datetime,
) -> None:
    """Plot all baseline lag peaks and their correlation coefficients."""
    time_values = [
        dt.datetime.fromtimestamp(value, tz=dt.timezone.utc)
        for value in wind_times
    ]
    channel_names = ("ch1", "ch3", "ch4")
    figure, axes = plt.subplots(
        len(fading_wind.PAIR_INDICES),
        2,
        figsize=(16, 13),
        sharex=True,
        sharey=True,
        facecolor="#070b14",
        constrained_layout=True,
    )
    for pair_index, (first, second) in enumerate(
        fading_wind.PAIR_INDICES
    ):
        baseline = (
            fading_wind.ANTENNA_EN_M[second]
            - fading_wind.ANTENNA_EN_M[first]
        )
        baseline_length = float(np.linalg.norm(baseline))
        label = (
            f"{channel_names[first]} → {channel_names[second]} · "
            f"ΔE={baseline[0]:+.1f} m, ΔN={baseline[1]:+.1f} m · "
            f"{baseline_length:.1f} m"
        )
        delay_image = axes[pair_index, 0].pcolormesh(
            time_values,
            fading_wind.ALTITUDE_KM,
            baseline_delay_s[:, :, pair_index].T,
            shading="nearest",
            cmap="seismic",
            vmin=-fading_wind.MAX_LAG_S,
            vmax=fading_wind.MAX_LAG_S,
        )
        axes[pair_index, 0].set_title(
            f"{label}\npeak delay",
            color="#edf4ff",
        )
        delay_colorbar = figure.colorbar(
            delay_image,
            ax=axes[pair_index, 0],
            pad=0.01,
        )
        delay_colorbar.set_label(
            "Delay (s)",
            color="#b7c5d9",
        )
        delay_colorbar.ax.tick_params(colors="#8fa1ba")

        correlation_image = axes[pair_index, 1].pcolormesh(
            time_values,
            fading_wind.ALTITUDE_KM,
            baseline_peak_correlation[:, :, pair_index].T,
            shading="nearest",
            cmap="viridis",
            vmin=0.0,
            vmax=1.0,
        )
        axes[pair_index, 1].set_title(
            f"{label}\npeak correlation",
            color="#edf4ff",
        )
        correlation_colorbar = figure.colorbar(
            correlation_image,
            ax=axes[pair_index, 1],
            pad=0.01,
        )
        correlation_colorbar.set_label(
            "Correlation coefficient",
            color="#b7c5d9",
        )
        correlation_colorbar.ax.tick_params(colors="#8fa1ba")

    for axis in axes.ravel():
        axis.set_facecolor("#050810")
        axis.set_xlim(start, end)
        axis.set_ylim(
            fading_wind.ALTITUDE_KM[0]
            - fading_wind.ALTITUDE_HALF_WIDTH_KM,
            fading_wind.ALTITUDE_KM[-1]
            + fading_wind.ALTITUDE_HALF_WIDTH_KM,
        )
        axis.set_ylabel("Altitude (km)", color="#b7c5d9")
        axis.tick_params(colors="#8fa1ba")
        for spine in axis.spines.values():
            spine.set_color("#314563")
    for axis in axes[-1]:
        axis.set_xlabel("Time (UTC)", color="#b7c5d9")
        axis.xaxis.set_major_locator(mdates.MinuteLocator(interval=5))
        axis.xaxis.set_major_formatter(
            mdates.DateFormatter("%H:%M", tz=dt.timezone.utc)
        )
    figure.suptitle(
        "Ramfjordmoen MF radar · fading-correlation quality audit · "
        "latest 15 minutes",
        color="#edf4ff",
        fontsize=18,
        weight="semibold",
    )
    temporary = FADING_DIAGNOSTIC_PLOT.with_suffix(".tmp.png")
    figure.savefig(
        temporary,
        dpi=150,
        facecolor=figure.get_facecolor(),
    )
    plt.close(figure)
    os.replace(temporary, FADING_DIAGNOSTIC_PLOT)


def main() -> None:
    reader = drf.DigitalMetadataReader(mc.xc_dir)
    bounds = reader.get_bounds()
    end_unix = float(bounds[1]) / 1e6 + 2.0
    start_unix = end_unix - WINDOW_S
    start = dt.datetime.fromtimestamp(start_unix, tz=dt.timezone.utc)
    end = dt.datetime.fromtimestamp(end_unix, tz=dt.timezone.utc)

    channel_times, channel_ranges, channel_power = load_channel_power(
        reader,
        start_unix,
        end_unix,
    )
    noise_mask = (
        (channel_ranges >= monitor.THIRTY_MINUTE_NOISE_RANGE_KM[0])
        & (channel_ranges <= monitor.THIRTY_MINUTE_NOISE_RANGE_KM[1])
    )
    broadband_noise_power = np.nanmean(
        channel_power[:, :, noise_mask],
        axis=(0, 2),
    )
    doppler_noise_power = np.maximum(
        broadband_noise_power[DOPPLER_CHANNEL_INDICES],
        1e-20,
    )
    (
        times,
        ranges,
        frequency_hz,
        velocity,
        fit_snr,
        fitted_channel_power,
    ) = fit_ten_second_doppler(
        reader,
        start_unix,
        end_unix,
    )
    narrowband_power_ratio = np.sum(
        fitted_channel_power / doppler_noise_power[None, :, None],
        axis=1,
    )
    if not np.array_equal(times, channel_times):
        raise RuntimeError("SNR and Doppler ten-second timestamps differ")
    if not np.array_equal(ranges, channel_ranges):
        raise RuntimeError("SNR and Doppler range vectors differ")
    (
        wind_times,
        zonal_wind_ms,
        meridional_wind_ms,
        wind_peak_correlation,
        wind_fit_rmse,
        baseline_delay_s,
        baseline_peak_correlation,
    ) = fading_wind.estimate_winds(
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
    plot_fading_correlation_diagnostics(
        wind_times,
        baseline_delay_s,
        baseline_peak_correlation,
        start,
        end,
    )
    plot_combined_rti(
        times,
        ranges,
        narrowband_power_ratio,
        times,
        ranges,
        velocity,
        wind_times,
        zonal_wind_ms,
        meridional_wind_ms,
        start,
        end,
    )
    fading_history.update_history(
        reader,
        end_unix,
        max_windows=4,
    )
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).with_name("plot_fading_ccf_diagnostics.py")),
        ],
        check=True,
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
        handle.attrs["measurement_selection"] = (
            "none_every_0_to_300_km_time_range_cell_is_fitted"
        )
        handle.attrs["narrowband_power"] = (
            "sum_of_squared_independent_complex_dipole_amplitudes"
        )
        handle.attrs["broadband_noise_power"] = (
            "mean_processed_voltage_power_over_30_to_50_km_and_15_minutes"
        )
        handle.attrs["fading_wind_method"] = (
            "evolving_elliptical_full_correlation_analysis"
        )
        handle.attrs["fading_wind_pattern_to_neutral_factor"] = 0.5
        handle.attrs["fading_wind_window_seconds"] = fading_wind.WINDOW_S
        handle.attrs["fading_wind_output_cadence_seconds"] = (
            fading_wind.OUTPUT_CADENCE_S
        )
        handle.attrs["fading_baseline_pairs"] = (
            "ch1_to_ch3,ch1_to_ch4,ch3_to_ch4"
        )
        handle.attrs["fading_delay_sign"] = (
            "positive_when_second_named_antenna_sees_pattern_later"
        )
        handle.attrs["fading_delay_search_limit_seconds"] = (
            fading_wind.MAX_LAG_S
        )
        handle.create_dataset("time_unix", data=times)
        handle.create_dataset("range_km", data=ranges)
        handle.create_dataset(
            "frequency_hz",
            data=frequency_hz,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
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
        handle.create_dataset(
            "processed_broadband_noise_power_by_channel",
            data=doppler_noise_power,
        )
        handle.create_dataset(
            "narrowband_to_broadband_noise_ratio",
            data=narrowband_power_ratio,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
        handle.create_dataset("wind_time_unix", data=wind_times)
        handle.create_dataset(
            "wind_altitude_km",
            data=fading_wind.ALTITUDE_KM,
        )
        handle.create_dataset(
            "zonal_wind_ms",
            data=zonal_wind_ms,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
        handle.create_dataset(
            "meridional_wind_ms",
            data=meridional_wind_ms,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
        handle.create_dataset(
            "wind_peak_correlation",
            data=wind_peak_correlation,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
        handle.create_dataset(
            "wind_fit_rmse",
            data=wind_fit_rmse,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
        handle.create_dataset(
            "fading_baseline_delay_s",
            data=baseline_delay_s,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
        handle.create_dataset(
            "fading_baseline_peak_correlation",
            data=baseline_peak_correlation,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
    os.replace(temporary_data, DATA_FILE)
    print(f"Four-channel SNR health: {CHANNEL_HEALTH_PLOT}")
    print(f"Combined SNR/Doppler RTI: {COMBINED_PLOT}")
    print(f"Fading correlation diagnostics: {FADING_DIAGNOSTIC_PLOT}")
    print(f"Dense fit cells: {velocity.size}")


if __name__ == "__main__":
    main()
