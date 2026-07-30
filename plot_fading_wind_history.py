#!/usr/bin/env python3
"""Maintain and plot a rolling 48-hour fading-wind history."""

from __future__ import annotations

import argparse
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

import fading_wind
import mf_conf as mc


HISTORY_HOURS = 48
WINDOW_S = 15 * 60
CADENCE_S = 15 * 60
PROCESSING_VERSION = 1
DISPLAY_LIMIT_MS = 150.0
DEFAULT_STATE = Path("/data2/plots/monitor/fading_wind_48h_state.h5")
DEFAULT_OUTPUT_DIR = Path("/data2/plots/monitor")
WIND_PLOT_NAME = "latest_fading_wind_48h.png"
QUALITY_PLOT_NAME = "latest_fading_wind_quality_48h.png"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-dir", default=mc.xc_dir)
    parser.add_argument("--state-file", type=Path, default=DEFAULT_STATE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--max-windows",
        type=int,
        default=4,
        help="Missing 15-minute windows to process; 0 processes all.",
    )
    return parser.parse_args()


def empty_state(times: np.ndarray) -> dict[str, np.ndarray]:
    shape = (len(times), len(fading_wind.ALTITUDE_KM))
    baseline_shape = (*shape, len(fading_wind.PAIR_INDICES))
    return {
        "time_unix": times,
        "zonal_wind_ms": np.full(shape, np.nan),
        "meridional_wind_ms": np.full(shape, np.nan),
        "wind_peak_correlation": np.full(shape, np.nan),
        "wind_fit_rmse": np.full(shape, np.nan),
        "baseline_delay_s": np.full(baseline_shape, np.nan),
        "baseline_peak_correlation": np.full(baseline_shape, np.nan),
    }


def load_state(path: Path, desired_times: np.ndarray) -> dict[str, np.ndarray]:
    result = empty_state(desired_times)
    if not path.exists():
        return result
    try:
        with h5py.File(path, "r") as handle:
            if int(handle.attrs["processing_version"]) != PROCESSING_VERSION:
                raise ValueError("processing version changed")
            old_times = np.asarray(handle["time_unix"], dtype=np.int64)
            old_altitude = np.asarray(handle["altitude_km"], dtype=np.float64)
            if not np.array_equal(old_altitude, fading_wind.ALTITUDE_KM):
                raise ValueError("altitude grid changed")
            shared, new_index, old_index = np.intersect1d(
                desired_times,
                old_times,
                return_indices=True,
            )
            del shared
            for name in result:
                if name == "time_unix":
                    continue
                old_value = np.asarray(handle[name], dtype=np.float64)
                result[name][new_index] = old_value[old_index]
        return result
    except Exception as error:
        print(f"Ignoring unreadable fading-wind state {path}: {error}")
        return empty_state(desired_times)


def atomic_state(path: Path, state: dict[str, np.ndarray]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["processing_version"] = PROCESSING_VERSION
        handle.attrs["integration_seconds"] = WINDOW_S
        handle.attrs["output_cadence_seconds"] = CADENCE_S
        handle.attrs["correlation_sample_rate_hz"] = (
            fading_wind.CORRELATION_SAMPLE_RATE_HZ
        )
        handle.attrs["pattern_to_neutral_wind_factor"] = 0.5
        handle.attrs["baseline_pairs"] = (
            "ch1_to_ch3,ch1_to_ch4,ch3_to_ch4"
        )
        handle.attrs["delay_sign"] = (
            "positive_when_second_named_antenna_sees_pattern_later"
        )
        handle.create_dataset("time_unix", data=state["time_unix"])
        handle.create_dataset(
            "altitude_km",
            data=fading_wind.ALTITUDE_KM,
        )
        for name, value in state.items():
            if name == "time_unix":
                continue
            handle.create_dataset(
                name,
                data=value,
                compression="gzip",
                compression_opts=1,
                shuffle=True,
            )
    os.replace(temporary, path)


def estimate_window(reader, center_unix: int) -> dict[str, np.ndarray]:
    half_window = WINDOW_S / 2.0
    sample_times, power = fading_wind.load_altitude_power(
        reader,
        center_unix - half_window,
        center_unix + half_window,
    )
    order = np.argsort(sample_times)
    sample_times = sample_times[order]
    power = power[order]
    sample_interval_s = float(np.median(np.diff(sample_times)))
    expected_samples = int(round(WINDOW_S / sample_interval_s))
    if len(sample_times) < 0.9 * expected_samples:
        raise RuntimeError("less than 90 percent of fading samples available")

    result = {
        "zonal_wind_ms": np.full(len(fading_wind.ALTITUDE_KM), np.nan),
        "meridional_wind_ms": np.full(len(fading_wind.ALTITUDE_KM), np.nan),
        "wind_peak_correlation": np.full(
            len(fading_wind.ALTITUDE_KM),
            np.nan,
        ),
        "wind_fit_rmse": np.full(len(fading_wind.ALTITUDE_KM), np.nan),
        "baseline_delay_s": np.full(
            (len(fading_wind.ALTITUDE_KM), len(fading_wind.PAIR_INDICES)),
            np.nan,
        ),
        "baseline_peak_correlation": np.full(
            (len(fading_wind.ALTITUDE_KM), len(fading_wind.PAIR_INDICES)),
            np.nan,
        ),
    }
    for altitude_index in range(len(fading_wind.ALTITUDE_KM)):
        values = power[:, :, altitude_index]
        if not np.all(np.isfinite(values)):
            continue
        (
            result["baseline_delay_s"][altitude_index],
            result["baseline_peak_correlation"][altitude_index],
        ) = fading_wind.cross_correlation_delays(
            values,
            sample_interval_s,
        )
        (
            result["zonal_wind_ms"][altitude_index],
            result["meridional_wind_ms"][altitude_index],
            result["wind_peak_correlation"][altitude_index],
            result["wind_fit_rmse"][altitude_index],
        ) = fading_wind.fit_fading_window(
            values,
            sample_interval_s,
        )
    return result


def style_axis(axis) -> None:
    axis.set_facecolor("#050810")
    axis.tick_params(colors="#8fa1ba")
    axis.yaxis.label.set_color("#b7c5d9")
    axis.title.set_color("#edf4ff")
    for spine in axis.spines.values():
        spine.set_color("#314563")


def configure_time_axis(axis) -> None:
    axis.xaxis.set_major_locator(mdates.HourLocator(interval=6))
    axis.xaxis.set_major_formatter(
        mdates.DateFormatter("%d %b\n%H:%M", tz=dt.timezone.utc)
    )


def plot_winds(state: dict[str, np.ndarray], output_path: Path) -> None:
    times = [
        dt.datetime.fromtimestamp(value, tz=dt.timezone.utc)
        for value in state["time_unix"]
    ]
    figure, axes = plt.subplots(
        2,
        1,
        figsize=(16, 10),
        sharex=True,
        sharey=True,
        facecolor="#070b14",
        constrained_layout=True,
    )
    panels = (
        ("zonal_wind_ms", "Zonal neutral wind · east positive"),
        ("meridional_wind_ms", "Meridional neutral wind · north positive"),
    )
    for axis, (name, title) in zip(axes, panels):
        style_axis(axis)
        image = axis.pcolormesh(
            times,
            fading_wind.ALTITUDE_KM,
            state[name].T,
            shading="nearest",
            cmap="seismic",
            vmin=-DISPLAY_LIMIT_MS,
            vmax=DISPLAY_LIMIT_MS,
        )
        axis.set_title(title)
        axis.set_ylabel("Altitude (km)")
        colorbar = figure.colorbar(image, ax=axis, pad=0.01)
        colorbar.set_label("Wind velocity (m/s)", color="#b7c5d9")
        colorbar.ax.tick_params(colors="#8fa1ba")
    axes[-1].set_xlabel("Window center (UTC)", color="#b7c5d9")
    configure_time_axis(axes[-1])
    figure.suptitle(
        "Ramfjordmoen MF radar · 15-minute fading-wind estimates · "
        "latest 48 hours",
        color="#edf4ff",
        fontsize=18,
        weight="semibold",
    )
    temporary = output_path.with_suffix(".tmp.png")
    figure.savefig(temporary, dpi=150, facecolor=figure.get_facecolor())
    plt.close(figure)
    os.replace(temporary, output_path)


def plot_quality(state: dict[str, np.ndarray], output_path: Path) -> None:
    times = [
        dt.datetime.fromtimestamp(value, tz=dt.timezone.utc)
        for value in state["time_unix"]
    ]
    channel_names = ("ch1", "ch3", "ch4")
    figure, axes = plt.subplots(
        4,
        2,
        figsize=(17, 17),
        sharex=True,
        sharey=True,
        facecolor="#070b14",
        constrained_layout=True,
    )
    for pair_index, (first, second) in enumerate(
        fading_wind.PAIR_INDICES
    ):
        label = f"{channel_names[first]} → {channel_names[second]}"
        delay_image = axes[pair_index, 0].pcolormesh(
            times,
            fading_wind.ALTITUDE_KM,
            state["baseline_delay_s"][:, :, pair_index].T,
            shading="nearest",
            cmap="seismic",
            vmin=-fading_wind.MAX_LAG_S,
            vmax=fading_wind.MAX_LAG_S,
        )
        axes[pair_index, 0].set_title(f"{label} peak delay")
        delay_colorbar = figure.colorbar(
            delay_image,
            ax=axes[pair_index, 0],
            pad=0.01,
        )
        delay_colorbar.set_label("Delay (s)", color="#b7c5d9")
        delay_colorbar.ax.tick_params(colors="#8fa1ba")

        correlation_image = axes[pair_index, 1].pcolormesh(
            times,
            fading_wind.ALTITUDE_KM,
            state["baseline_peak_correlation"][:, :, pair_index].T,
            shading="nearest",
            cmap="viridis",
            vmin=0.0,
            vmax=1.0,
        )
        axes[pair_index, 1].set_title(f"{label} peak correlation")
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

    rmse_image = axes[3, 0].pcolormesh(
        times,
        fading_wind.ALTITUDE_KM,
        state["wind_fit_rmse"].T,
        shading="nearest",
        cmap="magma",
        vmin=0.0,
        vmax=0.5,
    )
    axes[3, 0].set_title("Evolving-ellipse normalized fit RMSE")
    rmse_colorbar = figure.colorbar(rmse_image, ax=axes[3, 0], pad=0.01)
    rmse_colorbar.set_label("RMSE", color="#b7c5d9")
    rmse_colorbar.ax.tick_params(colors="#8fa1ba")

    accepted = np.isfinite(state["zonal_wind_ms"]).astype(np.float64)
    accepted_image = axes[3, 1].pcolormesh(
        times,
        fading_wind.ALTITUDE_KM,
        accepted.T,
        shading="nearest",
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
    )
    axes[3, 1].set_title("Automatic wind acceptance")
    accepted_colorbar = figure.colorbar(
        accepted_image,
        ax=axes[3, 1],
        pad=0.01,
        ticks=(0.0, 1.0),
    )
    accepted_colorbar.ax.set_yticklabels(("rejected", "accepted"))
    accepted_colorbar.ax.tick_params(colors="#8fa1ba")

    for axis in axes.ravel():
        style_axis(axis)
        axis.set_ylabel("Altitude (km)")
    for axis in axes[-1]:
        axis.set_xlabel("Window center (UTC)", color="#b7c5d9")
        configure_time_axis(axis)
    figure.suptitle(
        "Ramfjordmoen MF radar · 15-minute fading-wind quality panels · "
        "latest 48 hours",
        color="#edf4ff",
        fontsize=18,
        weight="semibold",
    )
    temporary = output_path.with_suffix(".tmp.png")
    figure.savefig(temporary, dpi=150, facecolor=figure.get_facecolor())
    plt.close(figure)
    os.replace(temporary, output_path)


def update_history(
    reader,
    end_unix: float,
    state_path: Path = DEFAULT_STATE,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    max_windows: int = 4,
) -> dict[str, np.ndarray]:
    half_window = WINDOW_S // 2
    start_unix = end_unix - HISTORY_HOURS * 3600
    first_center = int(
        np.ceil((start_unix + half_window) / CADENCE_S) * CADENCE_S
    )
    final_center = int(
        np.floor((end_unix - half_window) / CADENCE_S) * CADENCE_S
    )
    desired_times = np.arange(
        first_center,
        final_center + 1,
        CADENCE_S,
        dtype=np.int64,
    )
    state = load_state(state_path, desired_times)
    completed = np.any(
        np.isfinite(state["baseline_peak_correlation"]),
        axis=(1, 2),
    )
    missing_indices = np.flatnonzero(~completed)[::-1]
    if max_windows > 0:
        missing_indices = missing_indices[:max_windows]
    added = 0
    unavailable = 0
    for index in missing_indices:
        center = int(desired_times[index])
        try:
            result = estimate_window(reader, center)
        except Exception:
            unavailable += 1
            continue
        for name, value in result.items():
            state[name][index] = value
        added += 1

    state_path.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    atomic_state(state_path, state)
    plot_winds(state, output_dir / WIND_PLOT_NAME)
    plot_quality(state, output_dir / QUALITY_PLOT_NAME)
    print(
        f"Fading-wind history: calculated {added}, "
        f"unavailable {unavailable}, "
        f"remaining {np.count_nonzero(~np.any(np.isfinite(state['baseline_peak_correlation']), axis=(1, 2)))}"
    )
    return state


def main() -> None:
    args = parse_args()
    reader = drf.DigitalMetadataReader(args.metadata_dir)
    end_unix = float(reader.get_bounds()[1]) / 1e6
    update_history(
        reader,
        end_unix,
        args.state_file,
        args.output_dir,
        args.max_windows,
    )


if __name__ == "__main__":
    main()
