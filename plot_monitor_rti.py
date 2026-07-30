#!/usr/bin/env python3
"""Generate rolling 48-hour RTI quick-look plots for the web monitor."""

from __future__ import annotations

import argparse
import datetime as dt
import h5py
import json
import os
from pathlib import Path
import time

import digital_rf as drf
import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as constants

import mf_conf as mc
import radar_common as common
from plot_dense_doppler import fit_common_sinusoid_fft


FS = 1_000_000
IPP_SAMPLES = 10_000
DECIMATION = 10
OFFSET = 8_900
TX_REFERENCE_CHANNEL = "ch1"
TX_REFERENCE_SAMPLES = 100
GROUND_CLUTTER_SAMPLES = 120
TX_CENTER_RAW_SAMPLES = 54
FIR_TAPS = 50
DEFAULT_CHANNELS = ("ch1", "ch3", "ch4")
NOISE_RANGE_KM = (250.0, 1400.0)
NOISE_PERCENTILE = 20.0
THIRTY_MINUTE_NOISE_RANGE_KM = (30.0, 50.0)
THIRTY_MINUTE_POWER_RATIO_MIN_DB = -10.0
THIRTY_MINUTE_CADENCE_S = 2
THIRTY_MINUTE_WINDOW_S = 30 * 60
PROCESSING_VERSION = 2
DOPPLER_PROCESSING_VERSION = 1
DOPPLER_CADENCE_S = 60
DOPPLER_FIT_DURATION_S = 10
DOPPLER_SOURCE_RECORD_S = 2
DOPPLER_RANGE_MAX_KM = 200.0
DOPPLER_DISPLAY_LIMIT_MS = 150.0
DOPPLER_FIELDS = ("rti1", "rti3", "rti4")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", default=mc.raw_dir)
    parser.add_argument("--output-dir", default="/data2/plots/monitor")
    parser.add_argument("--state-file", default="/data2/plots/monitor/rti_48h_state.h5")
    parser.add_argument(
        "--doppler-state-file",
        default="/data2/plots/monitor/doppler_48h_state.h5",
    )
    parser.add_argument(
        "--thirty-minute-state-file",
        default="/data2/plots/monitor/rti_30m_2s_state.h5",
    )
    parser.add_argument("--hours", type=float, default=48.0)
    parser.add_argument("--cadence", type=int, default=60, help="Seconds per time bin.")
    parser.add_argument(
        "--averages",
        type=int,
        default=5,
        help="Consecutive IPPs averaged into each time bin.",
    )
    parser.add_argument("--channels", nargs="+", default=list(DEFAULT_CHANNELS))
    return parser.parse_args()


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def atomic_state(path: Path, times: np.ndarray, power: np.ndarray) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["processing_version"] = PROCESSING_VERSION
        handle.create_dataset(
            "times",
            data=times,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
        handle.create_dataset(
            "power",
            data=power,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
    os.replace(temporary, path)


def load_state(path: Path, n_range: int) -> tuple[np.ndarray, np.ndarray]:
    if not path.exists():
        return np.empty(0, dtype=np.int64), np.empty((0, n_range), dtype=np.float32)
    try:
        with h5py.File(path, "r") as state:
            version = int(state.attrs["processing_version"])
            if version != PROCESSING_VERSION:
                raise ValueError(
                    f"processing version {version} is not {PROCESSING_VERSION}"
                )
            times = np.asarray(state["times"], dtype=np.int64)
            power = np.asarray(state["power"], dtype=np.float32)
        if power.ndim != 2 or power.shape[1] != n_range or len(times) != len(power):
            raise ValueError("state dimensions do not match current RTI settings")
        return times, power
    except Exception as error:
        print(f"Ignoring unreadable RTI state {path}: {error}")
        return np.empty(0, dtype=np.int64), np.empty((0, n_range), dtype=np.float32)


def atomic_doppler_state(
    path: Path,
    times: np.ndarray,
    range_km: np.ndarray,
    velocity_ms: np.ndarray,
) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["processing_version"] = DOPPLER_PROCESSING_VERSION
        handle.attrs["fit_duration_seconds"] = DOPPLER_FIT_DURATION_S
        handle.attrs["output_cadence_seconds"] = DOPPLER_CADENCE_S
        handle.attrs["measurement_selection"] = "none"
        handle.create_dataset("times", data=times)
        handle.create_dataset("range_km", data=range_km)
        handle.create_dataset(
            "velocity_ms",
            data=velocity_ms,
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
    os.replace(temporary, path)


def load_doppler_state(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    empty = (
        np.empty(0, dtype=np.int64),
        np.empty(0, dtype=np.float64),
        np.empty((0, 0), dtype=np.float32),
    )
    if not path.exists():
        return empty
    try:
        with h5py.File(path, "r") as state:
            version = int(state.attrs["processing_version"])
            if version != DOPPLER_PROCESSING_VERSION:
                raise ValueError(
                    f"processing version {version} is not "
                    f"{DOPPLER_PROCESSING_VERSION}"
                )
            times = np.asarray(state["times"], dtype=np.int64)
            range_km = np.asarray(state["range_km"], dtype=np.float64)
            velocity_ms = np.asarray(state["velocity_ms"], dtype=np.float32)
        if (
            velocity_ms.ndim != 2
            or velocity_ms.shape != (len(times), len(range_km))
        ):
            raise ValueError("Doppler state dimensions are inconsistent")
        return times, range_km, velocity_ms
    except Exception as error:
        print(f"Ignoring unreadable Doppler state {path}: {error}")
        return empty


def process_doppler_bin(
    reader: drf.DigitalMetadataReader,
    group_start: int,
) -> tuple[np.ndarray, np.ndarray]:
    expected_keys = [
        int((group_start + offset) * 1e6)
        for offset in range(
            0,
            DOPPLER_FIT_DURATION_S,
            DOPPLER_SOURCE_RECORD_S,
        )
    ]
    metadata = reader.read(
        expected_keys[0],
        int((group_start + DOPPLER_FIT_DURATION_S) * 1e6) - 1,
    )
    if not all(key in metadata for key in expected_keys):
        raise RuntimeError("incomplete ten-second metadata group")
    records = [metadata[key] for key in expected_keys]
    required = (*DOPPLER_FIELDS, "rvec", "tvec")
    if not all(field in record for record in records for field in required):
        raise RuntimeError("Doppler metadata fields are incomplete")

    range_km = np.asarray(records[0]["rvec"], dtype=np.float64)
    if not all(
        np.array_equal(
            range_km,
            np.asarray(record["rvec"], dtype=np.float64),
        )
        for record in records[1:]
    ):
        raise RuntimeError("range vectors differ within Doppler group")
    mask = (range_km >= 0.0) & (range_km <= DOPPLER_RANGE_MAX_KM)
    selected_range = range_km[mask]
    fit_times = np.concatenate(
        [
            float(key) / 1e6
            + np.asarray(record["tvec"], dtype=np.float64)
            for key, record in zip(expected_keys, records)
        ]
    )
    fit_times -= fit_times[0]
    channel_voltage = np.asarray(
        [
            np.concatenate(
                [
                    np.asarray(record[field])[:, mask]
                    for record in records
                ],
                axis=0,
            )
            for field in DOPPLER_FIELDS
        ],
        dtype=np.complex128,
    )
    frequency_hz, _, _ = fit_common_sinusoid_fft(
        fit_times,
        channel_voltage.transpose(1, 0, 2),
        common.MAX_FIT_DOPPLER_HZ,
    )
    velocity_ms = (
        common.DOPPLER_SIGN
        * common.DOPPLER_TO_MS
        * frequency_hz
    )
    return selected_range, velocity_ms.astype(np.float32)


def process_time_bin(
    reader: drf.DigitalRFReader,
    unix_time: int,
    channels: list[str],
    averages: int,
) -> np.ndarray:
    lpf = mc.fir_lowpass_hann(fc=20e3, fs=FS, num_taps=FIR_TAPS)
    decimated_indices = np.arange(IPP_SAMPLES // DECIMATION) * DECIMATION
    coherent_voltage = {
        channel: np.zeros(IPP_SAMPLES // DECIMATION, dtype=np.complex64)
        for channel in channels
    }
    first_sample = unix_time * FS + OFFSET
    block_samples = averages * IPP_SAMPLES
    voltage_blocks = {
        channel: (
            reader.read_vector_1d(first_sample, block_samples, channel).astype(
                np.complex64,
                casting="unsafe",
                copy=False,
            )
            - mc.dc_offset
        ).reshape(averages, IPP_SAMPLES)
        for channel in channels
    }
    tx_block = (
        voltage_blocks[TX_REFERENCE_CHANNEL]
        if TX_REFERENCE_CHANNEL in voltage_blocks
        else (
            reader.read_vector_1d(
                first_sample,
                block_samples,
                TX_REFERENCE_CHANNEL,
            ).astype(np.complex64, casting="unsafe", copy=False)
            - mc.dc_offset
        ).reshape(averages, IPP_SAMPLES)
    )

    for ipp_index in range(averages):
        tx_voltage = tx_block[ipp_index]
        tx_phase = np.angle(np.mean(tx_voltage[:TX_REFERENCE_SAMPLES]))

        for channel in channels:
            voltage = voltage_blocks[channel][ipp_index].copy()
            voltage[:GROUND_CLUTTER_SAMPLES] = 0.0
            filtered = np.convolve(
                np.exp(-1j * tx_phase) * voltage,
                lpf,
                mode="same",
            )
            coherent_voltage[channel] += filtered[decimated_indices]

    channel_powers = [
        np.abs(coherent_voltage[channel] / averages) ** 2
        for channel in channels
    ]
    return np.mean(channel_powers, axis=0).astype(np.float32)


def range_vector_km(n_range: int) -> np.ndarray:
    raw_indices = np.arange(n_range, dtype=np.float64) * DECIMATION
    return (
        (raw_indices - TX_CENTER_RAW_SAMPLES)
        * constants.c
        / (2.0 * FS * 1_000.0)
    )


def plot_rti(
    times: np.ndarray,
    power: np.ndarray,
    output_path: Path,
    maximum_range_km: float,
    title: str,
    display_min_db: float,
    display_max_db: float,
    fixed_noise_power: float | None = None,
) -> None:
    range_km = range_vector_km(power.shape[1])
    mask = range_km <= maximum_range_km
    noise_mask = (range_km >= NOISE_RANGE_KM[0]) & (range_km <= NOISE_RANGE_KM[1])
    if fixed_noise_power is None:
        noise_power = np.nanpercentile(power[:, noise_mask], NOISE_PERCENTILE, axis=1)
    else:
        noise_power = np.full(power.shape[0], fixed_noise_power, dtype=np.float64)
    noise_power = np.maximum(noise_power, 1e-20)
    power_ratio = np.maximum(
        power[:, mask] / noise_power[:, np.newaxis],
        1e-20,
    )
    power_ratio_db = (10.0 * np.log10(power_ratio)).T
    power_ratio_db = np.maximum(power_ratio_db, display_min_db)
    finite = power_ratio_db[np.isfinite(power_ratio_db)]
    if finite.size == 0:
        raise RuntimeError("RTI contains no finite samples")

    date_numbers = mdates.date2num(times.astype("datetime64[s]"))
    half_step_days = (
        0.5 * float(np.median(np.diff(date_numbers)))
        if len(date_numbers) > 1
        else 0.5 / (24.0 * 60.0)
    )
    range_step_km = (
        float(np.median(np.diff(range_km[mask])))
        if np.count_nonzero(mask) > 1
        else DECIMATION * constants.c / (2.0 * FS * 1_000.0)
    )
    figure, axis = plt.subplots(figsize=(15, 5.8), facecolor="#070b14")
    axis.set_facecolor("#050810")
    mesh = axis.imshow(
        power_ratio_db,
        extent=(
            date_numbers[0] - half_step_days,
            date_numbers[-1] + half_step_days,
            range_km[mask][0] - 0.5 * range_step_km,
            range_km[mask][-1] + 0.5 * range_step_km,
        ),
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        cmap="plasma",
        vmin=display_min_db,
        vmax=display_max_db,
    )
    axis.xaxis_date()
    axis.set_ylim(0, maximum_range_km)
    axis.set_title(title, color="#edf4ff", fontsize=17, pad=14, weight="semibold")
    axis.set_xlabel("Time (UTC)", color="#b7c5d9", labelpad=9)
    axis.set_ylabel("Round-trip range (km)", color="#b7c5d9", labelpad=9)
    duration_seconds = float(times[-1] - times[0])
    if duration_seconds <= 3600:
        axis.xaxis.set_major_locator(mdates.MinuteLocator(interval=5))
        axis.xaxis.set_major_formatter(
            mdates.DateFormatter("%H:%M", tz=dt.timezone.utc)
        )
    else:
        axis.xaxis.set_major_locator(mdates.HourLocator(interval=6))
        axis.xaxis.set_major_formatter(
            mdates.DateFormatter("%d %b\n%H:%M", tz=dt.timezone.utc)
        )
    axis.tick_params(colors="#8fa1ba", labelsize=10)
    for spine in axis.spines.values():
        spine.set_color("#314563")
    axis.grid(color="#ffffff", alpha=0.07, linewidth=0.6)

    colorbar = figure.colorbar(mesh, ax=axis, pad=0.015)
    colorbar.set_label("Power / background (dB)", color="#b7c5d9")
    colorbar.ax.tick_params(colors="#8fa1ba")
    colorbar.outline.set_edgecolor("#314563")
    figure.tight_layout()

    temporary = output_path.with_suffix(".tmp.png")
    figure.savefig(temporary, dpi=150, facecolor=figure.get_facecolor(), bbox_inches="tight")
    plt.close(figure)
    os.replace(temporary, output_path)


def plot_mesosphere_power_doppler(
    power_times: np.ndarray,
    power: np.ndarray,
    doppler_times: np.ndarray,
    doppler_range_km: np.ndarray,
    velocity_ms: np.ndarray,
    output_path: Path,
) -> None:
    """Plot the 48-hour 0–200 km power RTI with an ungated Doppler panel."""
    power_range_km = range_vector_km(power.shape[1])
    power_mask = (
        (power_range_km >= 0.0)
        & (power_range_km <= DOPPLER_RANGE_MAX_KM)
    )
    noise_mask = (
        (power_range_km >= NOISE_RANGE_KM[0])
        & (power_range_km <= NOISE_RANGE_KM[1])
    )
    noise_power = np.nanpercentile(
        power[:, noise_mask],
        NOISE_PERCENTILE,
        axis=1,
    )
    power_ratio_db = 10.0 * np.log10(
        np.maximum(
            power[:, power_mask]
            / np.maximum(noise_power[:, None], 1e-20),
            1e-20,
        )
    )

    power_dates = mdates.date2num(power_times.astype("datetime64[s]"))
    doppler_dates = mdates.date2num(
        doppler_times.astype("datetime64[s]")
    )
    power_half_step = (
        0.5 * float(np.median(np.diff(power_dates)))
        if len(power_dates) > 1
        else 0.5 / (24.0 * 60.0)
    )
    doppler_half_step = (
        0.5 * float(np.median(np.diff(doppler_dates)))
        if len(doppler_dates) > 1
        else 0.5 / (24.0 * 60.0)
    )
    power_range_step = float(np.median(np.diff(power_range_km[power_mask])))
    doppler_range_step = float(np.median(np.diff(doppler_range_km)))

    figure, axes = plt.subplots(
        2,
        1,
        figsize=(15, 10.5),
        sharex=True,
        sharey=True,
        facecolor="#070b14",
        constrained_layout=True,
    )
    for axis in axes:
        axis.set_facecolor("#050810")
        axis.set_ylim(0.0, DOPPLER_RANGE_MAX_KM)
        axis.set_ylabel("Round-trip range (km)", color="#b7c5d9")
        axis.tick_params(colors="#8fa1ba", labelsize=10)
        axis.grid(color="#ffffff", alpha=0.07, linewidth=0.6)
        for spine in axis.spines.values():
            spine.set_color("#314563")

    power_image = axes[0].imshow(
        power_ratio_db.T,
        extent=(
            power_dates[0] - power_half_step,
            power_dates[-1] + power_half_step,
            power_range_km[power_mask][0] - 0.5 * power_range_step,
            power_range_km[power_mask][-1] + 0.5 * power_range_step,
        ),
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        cmap="plasma",
        vmin=-3.0,
        vmax=20.0,
    )
    axes[0].set_title(
        "Noise-referenced power",
        color="#edf4ff",
        fontsize=14,
    )
    power_colorbar = figure.colorbar(power_image, ax=axes[0], pad=0.015)
    power_colorbar.set_label(
        "Power / background (dB)",
        color="#b7c5d9",
    )
    power_colorbar.ax.tick_params(colors="#8fa1ba")
    power_colorbar.outline.set_edgecolor("#314563")

    doppler_image = axes[1].imshow(
        velocity_ms.T,
        extent=(
            doppler_dates[0] - doppler_half_step,
            doppler_dates[-1] + doppler_half_step,
            doppler_range_km[0] - 0.5 * doppler_range_step,
            doppler_range_km[-1] + 0.5 * doppler_range_step,
        ),
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        cmap="seismic",
        vmin=-DOPPLER_DISPLAY_LIMIT_MS,
        vmax=DOPPLER_DISPLAY_LIMIT_MS,
    )
    axes[1].set_title(
        "Three-dipole common-frequency 10-second Doppler fit",
        color="#edf4ff",
        fontsize=14,
    )
    doppler_colorbar = figure.colorbar(
        doppler_image,
        ax=axes[1],
        pad=0.015,
    )
    doppler_colorbar.set_label(
        "Monostatic Doppler velocity (m/s)",
        color="#b7c5d9",
    )
    doppler_colorbar.ax.tick_params(colors="#8fa1ba")
    doppler_colorbar.outline.set_edgecolor("#314563")

    axes[1].set_xlabel("Time (UTC)", color="#b7c5d9", labelpad=9)
    axes[1].xaxis.set_major_locator(mdates.HourLocator(interval=6))
    axes[1].xaxis.set_major_formatter(
        mdates.DateFormatter("%d %b\n%H:%M", tz=dt.timezone.utc)
    )
    axes[1].set_xlim(
        power_dates[0] - power_half_step,
        power_dates[-1] + power_half_step,
    )
    figure.suptitle(
        "Ramfjordmoen MF radar · latest 48 hours · 0–200 km",
        color="#edf4ff",
        fontsize=18,
        weight="semibold",
    )

    temporary = output_path.with_suffix(".tmp.png")
    figure.savefig(
        temporary,
        dpi=150,
        facecolor=figure.get_facecolor(),
    )
    plt.close(figure)
    os.replace(temporary, output_path)


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    state_path = Path(args.state_file)
    state_path.parent.mkdir(parents=True, exist_ok=True)

    reader = drf.DigitalRFReader(args.raw_dir)
    bounds = [reader.get_bounds(channel) for channel in args.channels]
    earliest = max(bound[0] for bound in bounds)
    latest = min(bound[1] for bound in bounds) - args.averages * IPP_SAMPLES - OFFSET
    if latest <= earliest:
        raise RuntimeError("No common raw-data interval is available across the selected channels")

    end_time = int(latest // FS // args.cadence * args.cadence)
    start_time = end_time - int(args.hours * 3600)
    raw_start_time = int(np.ceil(earliest / FS / args.cadence) * args.cadence)
    start_time = max(start_time, raw_start_time)

    n_range = IPP_SAMPLES // DECIMATION
    times, power = load_state(state_path, n_range)
    keep = (times >= start_time) & (times <= end_time)
    times = times[keep]
    power = power[keep]

    next_time = start_time if len(times) == 0 else int(times[-1]) + args.cadence
    new_times = []
    new_power = []
    started = time.monotonic()
    for index, unix_time in enumerate(range(next_time, end_time + 1, args.cadence), start=1):
        try:
            value = process_time_bin(reader, unix_time, args.channels, args.averages)
        except Exception as error:
            print(f"Skipping {unix_time}: {error}")
            continue
        new_times.append(unix_time)
        new_power.append(value)
        if index % 100 == 0:
            elapsed = time.monotonic() - started
            print(f"Processed {index} new time bins in {elapsed:.1f} seconds")

    if new_times:
        times = np.concatenate((times, np.asarray(new_times, dtype=np.int64)))
        power = np.concatenate((power, np.asarray(new_power, dtype=np.float32)), axis=0)

    if len(times) < 2:
        raise RuntimeError("Fewer than two RTI time bins are available")

    keep = times >= end_time - int(args.hours * 3600)
    times = times[keep]
    power = power[keep]
    atomic_state(state_path, times, power)

    doppler_reader = drf.DigitalMetadataReader(mc.xc_dir)
    doppler_bounds = doppler_reader.get_bounds()
    doppler_first_start = int(
        np.ceil(
            (float(doppler_bounds[0]) / 1e6)
            / DOPPLER_CADENCE_S
        )
        * DOPPLER_CADENCE_S
    )
    doppler_last_start = int(
        np.floor(
            (
                float(doppler_bounds[1]) / 1e6
                - DOPPLER_FIT_DURATION_S
            )
            / DOPPLER_CADENCE_S
        )
        * DOPPLER_CADENCE_S
    )
    doppler_window_start = max(
        end_time - int(args.hours * 3600),
        doppler_first_start,
    )
    doppler_state_path = Path(args.doppler_state_file)
    doppler_state_path.parent.mkdir(parents=True, exist_ok=True)
    doppler_times, doppler_range_km, velocity_ms = load_doppler_state(
        doppler_state_path
    )
    doppler_keep = (
        (doppler_times >= doppler_window_start)
        & (
            doppler_times
            <= doppler_last_start + DOPPLER_FIT_DURATION_S // 2
        )
    )
    doppler_times = doppler_times[doppler_keep]
    velocity_ms = velocity_ms[doppler_keep]
    next_doppler_start = (
        doppler_window_start
        if len(doppler_times) == 0
        else int(doppler_times[-1])
        - DOPPLER_FIT_DURATION_S // 2
        + DOPPLER_CADENCE_S
    )
    new_doppler_times = []
    new_velocity = []
    for index, group_start in enumerate(
        range(
            next_doppler_start,
            doppler_last_start + 1,
            DOPPLER_CADENCE_S,
        ),
        start=1,
    ):
        try:
            record_range, record_velocity = process_doppler_bin(
                doppler_reader,
                group_start,
            )
        except Exception as error:
            print(f"Skipping Doppler bin {group_start}: {error}")
            continue
        if len(doppler_range_km) == 0:
            doppler_range_km = record_range
        elif not np.array_equal(doppler_range_km, record_range):
            raise RuntimeError("Doppler range vector changed")
        new_doppler_times.append(
            group_start + DOPPLER_FIT_DURATION_S // 2
        )
        new_velocity.append(record_velocity)
        if index % 100 == 0:
            elapsed = time.monotonic() - started
            print(
                f"Processed {index} new Doppler bins in "
                f"{elapsed:.1f} seconds"
            )
    if new_doppler_times:
        doppler_times = np.concatenate(
            (
                doppler_times,
                np.asarray(new_doppler_times, dtype=np.int64),
            )
        )
        new_velocity_array = np.asarray(new_velocity, dtype=np.float32)
        velocity_ms = (
            new_velocity_array
            if velocity_ms.size == 0
            else np.concatenate((velocity_ms, new_velocity_array), axis=0)
        )
    if len(doppler_times) < 2:
        raise RuntimeError("Fewer than two 48-hour Doppler bins are available")
    atomic_doppler_state(
        doppler_state_path,
        doppler_times,
        doppler_range_km,
        velocity_ms,
    )

    plot_rti(
        times,
        power,
        output_dir / "latest_rti_48h_full.png",
        1500.0,
        "Ramfjordmoen MF radar · latest 48 hours · 0–1500 km",
        -3.0,
        35.0,
    )
    plot_mesosphere_power_doppler(
        times,
        power,
        doppler_times,
        doppler_range_km,
        velocity_ms,
        output_dir / "latest_rti_48h_mesosphere.png",
    )

    high_resolution_state = Path(args.thirty_minute_state_file)
    high_resolution_state.parent.mkdir(parents=True, exist_ok=True)
    high_resolution_end = int(
        latest // FS // THIRTY_MINUTE_CADENCE_S
        * THIRTY_MINUTE_CADENCE_S
    )
    high_resolution_start = max(
        high_resolution_end - THIRTY_MINUTE_WINDOW_S,
        int(
            np.ceil(earliest / FS / THIRTY_MINUTE_CADENCE_S)
            * THIRTY_MINUTE_CADENCE_S
        ),
    )
    recent_times, recent_power = load_state(
        high_resolution_state,
        n_range,
    )
    recent_keep = (
        (recent_times >= high_resolution_start)
        & (recent_times <= high_resolution_end)
    )
    recent_times = recent_times[recent_keep]
    recent_power = recent_power[recent_keep]
    next_recent_time = (
        high_resolution_start
        if len(recent_times) == 0
        else int(recent_times[-1]) + THIRTY_MINUTE_CADENCE_S
    )
    new_recent_times = []
    new_recent_power = []
    for unix_time in range(
        next_recent_time,
        high_resolution_end + 1,
        THIRTY_MINUTE_CADENCE_S,
    ):
        try:
            value = process_time_bin(
                reader,
                unix_time,
                args.channels,
                args.averages,
            )
        except Exception as error:
            print(f"Skipping two-second RTI bin {unix_time}: {error}")
            continue
        new_recent_times.append(unix_time)
        new_recent_power.append(value)
    if new_recent_times:
        recent_times = np.concatenate(
            (recent_times, np.asarray(new_recent_times, dtype=np.int64))
        )
        recent_power = np.concatenate(
            (
                recent_power,
                np.asarray(new_recent_power, dtype=np.float32),
            ),
            axis=0,
        )
    atomic_state(high_resolution_state, recent_times, recent_power)

    if len(recent_times) >= 2:
        range_km = range_vector_km(recent_power.shape[1])
        background_mask = (
            (range_km >= THIRTY_MINUTE_NOISE_RANGE_KM[0])
            & (range_km <= THIRTY_MINUTE_NOISE_RANGE_KM[1])
        )
        thirty_minute_noise_power = float(
            np.nanmean(recent_power[:, background_mask])
        )
    else:
        thirty_minute_noise_power = None

    generated_at = dt.datetime.now(dt.timezone.utc)
    status = {
        "generated_at": generated_at.isoformat().replace("+00:00", "Z"),
        "window_start": dt.datetime.fromtimestamp(int(times[0]), tz=dt.timezone.utc)
        .isoformat()
        .replace("+00:00", "Z"),
        "window_end": dt.datetime.fromtimestamp(int(times[-1]), tz=dt.timezone.utc)
        .isoformat()
        .replace("+00:00", "Z"),
        "time_bins": int(len(times)),
        "cadence_seconds": args.cadence,
        "averaged_ipps": args.averages,
        "channels": args.channels,
        "processing_version": PROCESSING_VERSION,
        "integration": "transmit_phase_referenced_coherent_voltage",
        "fir_cutoff_hz": 20_000,
        "fir_taps": FIR_TAPS,
        "range_resolution_km": DECIMATION * constants.c / (2.0 * FS * 1_000.0),
        "tx_center_raw_samples": TX_CENTER_RAW_SAMPLES,
        "noise_range_km": list(NOISE_RANGE_KM),
        "noise_percentile": NOISE_PERCENTILE,
        "display_quantity": "10log10(power/background_power)",
        "full_range_power_ratio_limits_db": [-3.0, 35.0],
        "mesosphere_power_ratio_limits_db": [-3.0, 20.0],
        "mesosphere_doppler_time_bins": int(len(doppler_times)),
        "mesosphere_doppler_cadence_seconds": DOPPLER_CADENCE_S,
        "mesosphere_doppler_fit_seconds": DOPPLER_FIT_DURATION_S,
        "mesosphere_doppler_channels": list(DOPPLER_FIELDS),
        "mesosphere_doppler_velocity_limits_ms": [
            -DOPPLER_DISPLAY_LIMIT_MS,
            DOPPLER_DISPLAY_LIMIT_MS,
        ],
        "mesosphere_doppler_measurement_selection": "none",
        "mesosphere_30m_time_bins": int(len(recent_times)),
        "mesosphere_30m_cadence_seconds": THIRTY_MINUTE_CADENCE_S,
        "mesosphere_30m_noise_method": "mean_power_over_time_and_range",
        "mesosphere_30m_noise_range_km": list(THIRTY_MINUTE_NOISE_RANGE_KM),
        "mesosphere_30m_noise_power": thirty_minute_noise_power,
        "mesosphere_30m_power_ratio_limits_db": [
            THIRTY_MINUTE_POWER_RATIO_MIN_DB,
            20.0,
        ],
    }
    atomic_json(output_dir / "rti_status.json", status)
    print(
        f"Updated {len(times)} RTI bins through {status['window_end']} "
        f"in {time.monotonic() - started:.1f} seconds"
    )


if __name__ == "__main__":
    main()
