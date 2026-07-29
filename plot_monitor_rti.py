#!/usr/bin/env python3
"""Generate rolling 48-hour RTI quick-look plots for the web monitor."""

from __future__ import annotations

import argparse
import datetime as dt
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


FS = 1_000_000
IPP_SAMPLES = 10_000
DECIMATION = 10
OFFSET = 8_900
DEFAULT_CHANNELS = ("ch1", "ch3", "ch4")
NOISE_RANGE_KM = (250.0, 1400.0)
NOISE_PERCENTILE = 20.0
THIRTY_MINUTE_NOISE_RANGE_KM = (30.0, 50.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", default=mc.raw_dir)
    parser.add_argument("--output-dir", default="/data2/plots/monitor")
    parser.add_argument("--state-file", default="/data2/plots/monitor/rti_48h_state.npz")
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
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, times=times, power=power)
    os.replace(temporary, path)


def load_state(path: Path, n_range: int) -> tuple[np.ndarray, np.ndarray]:
    if not path.exists():
        return np.empty(0, dtype=np.int64), np.empty((0, n_range), dtype=np.float32)

    try:
        with np.load(path) as state:
            times = state["times"].astype(np.int64, copy=False)
            power = state["power"].astype(np.float32, copy=False)
        if power.ndim != 2 or power.shape[1] != n_range or len(times) != len(power):
            raise ValueError("state dimensions do not match current RTI settings")
        return times, power
    except Exception as error:
        print(f"Ignoring unreadable RTI state {path}: {error}")
        return np.empty(0, dtype=np.int64), np.empty((0, n_range), dtype=np.float32)


def process_time_bin(
    reader: drf.DigitalRFReader,
    unix_time: int,
    channels: list[str],
    averages: int,
) -> np.ndarray:
    channel_powers = []
    first_sample = unix_time * FS + OFFSET

    for channel in channels:
        ipp_powers = []
        for ipp_index in range(averages):
            sample = first_sample + ipp_index * IPP_SAMPLES
            voltage = reader.read_vector_c81d(sample, IPP_SAMPLES, channel) - mc.dc_offset
            power = np.abs(voltage) ** 2
            ipp_powers.append(power.reshape(-1, DECIMATION).mean(axis=1))
        channel_powers.append(np.mean(ipp_powers, axis=0))

    return np.mean(channel_powers, axis=0).astype(np.float32)


def plot_rti(
    times: np.ndarray,
    power: np.ndarray,
    output_path: Path,
    maximum_range_km: float,
    title: str,
    snr_min_db: float,
    snr_max_db: float,
    fixed_noise_power: float | None = None,
) -> None:
    range_km = (
        np.arange(power.shape[1], dtype=np.float64)
        * DECIMATION
        * constants.c
        / (2.0 * FS * 1_000.0)
    )
    mask = range_km <= maximum_range_km
    noise_mask = (range_km >= NOISE_RANGE_KM[0]) & (range_km <= NOISE_RANGE_KM[1])
    if fixed_noise_power is None:
        noise_power = np.nanpercentile(power[:, noise_mask], NOISE_PERCENTILE, axis=1)
    else:
        noise_power = np.full(power.shape[0], fixed_noise_power, dtype=np.float64)
    noise_power = np.maximum(noise_power, 1e-20)
    signal_power = power[:, mask] - noise_power[:, np.newaxis]
    snr_linear = np.maximum(signal_power / noise_power[:, np.newaxis], 1e-20)
    snr_db = (10.0 * np.log10(snr_linear)).T
    snr_db = np.maximum(snr_db, snr_min_db)
    finite = snr_db[np.isfinite(snr_db)]
    if finite.size == 0:
        raise RuntimeError("RTI contains no finite samples")

    dates = [dt.datetime.fromtimestamp(int(value), tz=dt.timezone.utc) for value in times]
    figure, axis = plt.subplots(figsize=(15, 5.8), facecolor="#070b14")
    axis.set_facecolor("#050810")
    mesh = axis.pcolormesh(
        dates,
        range_km[mask],
        snr_db,
        cmap="plasma",
        shading="auto",
        vmin=snr_min_db,
        vmax=snr_max_db,
        rasterized=True,
    )
    axis.set_ylim(0, maximum_range_km)
    axis.set_title(title, color="#edf4ff", fontsize=17, pad=14, weight="semibold")
    axis.set_xlabel("Time (UTC)", color="#b7c5d9", labelpad=9)
    axis.set_ylabel("One-way range (km)", color="#b7c5d9", labelpad=9)
    axis.xaxis.set_major_locator(mdates.HourLocator(interval=6))
    axis.xaxis.set_major_formatter(mdates.DateFormatter("%d %b\n%H:%M", tz=dt.timezone.utc))
    axis.tick_params(colors="#8fa1ba", labelsize=10)
    for spine in axis.spines.values():
        spine.set_color("#314563")
    axis.grid(color="#ffffff", alpha=0.07, linewidth=0.6)

    colorbar = figure.colorbar(mesh, ax=axis, pad=0.015)
    colorbar.set_label("Quick-look SNR (dB)", color="#b7c5d9")
    colorbar.ax.tick_params(colors="#8fa1ba")
    colorbar.outline.set_edgecolor("#314563")
    figure.tight_layout()

    temporary = output_path.with_suffix(".tmp.png")
    figure.savefig(temporary, dpi=150, facecolor=figure.get_facecolor(), bbox_inches="tight")
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

    plot_rti(
        times,
        power,
        output_dir / "latest_rti_48h_full.png",
        1500.0,
        "Ramfjordmoen MF radar · latest 48 hours · 0–1500 km",
        -3.0,
        35.0,
    )
    plot_rti(
        times,
        power,
        output_dir / "latest_rti_48h_mesosphere.png",
        200.0,
        "Ramfjordmoen MF radar · latest 48 hours · 0–200 km",
        -3.0,
        20.0,
    )
    recent = times >= times[-1] - 30 * 60
    if np.count_nonzero(recent) >= 2:
        range_km = (
            np.arange(power.shape[1], dtype=np.float64)
            * DECIMATION
            * constants.c
            / (2.0 * FS * 1_000.0)
        )
        background_mask = (
            (range_km >= THIRTY_MINUTE_NOISE_RANGE_KM[0])
            & (range_km <= THIRTY_MINUTE_NOISE_RANGE_KM[1])
        )
        thirty_minute_noise_power = float(
            np.nanmean(power[recent][:, background_mask])
        )
        plot_rti(
            times[recent],
            power[recent],
            output_dir / "latest_rti_30m_mesosphere.png",
            200.0,
            "Ramfjordmoen MF radar · latest 30 minutes · 0–200 km",
            0.0,
            20.0,
            fixed_noise_power=thirty_minute_noise_power,
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
        "range_resolution_km": DECIMATION * constants.c / (2.0 * FS * 1_000.0),
        "noise_range_km": list(NOISE_RANGE_KM),
        "noise_percentile": NOISE_PERCENTILE,
        "full_range_snr_limits_db": [-3.0, 35.0],
        "mesosphere_snr_limits_db": [-3.0, 20.0],
        "mesosphere_30m_time_bins": int(np.count_nonzero(recent)),
        "mesosphere_30m_noise_method": "mean_power_over_time_and_range",
        "mesosphere_30m_noise_range_km": list(THIRTY_MINUTE_NOISE_RANGE_KM),
        "mesosphere_30m_noise_power": thirty_minute_noise_power,
        "mesosphere_30m_snr_limits_db": [0.0, 20.0],
    }
    atomic_json(output_dir / "rti_status.json", status)
    print(
        f"Updated {len(times)} RTI bins through {status['window_end']} "
        f"in {time.monotonic() - started:.1f} seconds"
    )


if __name__ == "__main__":
    main()
