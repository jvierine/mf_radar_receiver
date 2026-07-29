#!/usr/bin/env python3
"""Analyze a UTC interval with coherent RTIs and automatic wind detections."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path

import digital_rf as drf
import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

import mf_conf as mc
import plot_monitor_rti as monitor
import wind_estimates_3ch_2days as winds


DEFAULT_OUTPUT_ROOT = "/data2/products/interval_analysis"
WIND_BLOCK_S = 10 * 60
WIND_SEED_S = 10 * 60
POWER_THRESHOLD_DB = 6.0


def parse_utc(value: str) -> dt.datetime:
    parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    return parsed.astimezone(dt.timezone.utc)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("start", type=parse_utc, help="UTC interval start.")
    parser.add_argument("end", type=parse_utc, help="UTC interval end.")
    parser.add_argument("--raw-dir", default=mc.raw_dir)
    parser.add_argument("--metadata-dir", default=mc.xc_dir)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--cadence", type=int, default=60)
    parser.add_argument("--averages", type=int, default=5)
    return parser.parse_args()


def band_statistics(
    power_ratio_db: np.ndarray,
    range_km: np.ndarray,
    lower_km: float,
    upper_km: float,
) -> dict:
    mask = (range_km >= lower_km) & (range_km < upper_km)
    values = power_ratio_db[:, mask]
    profile_peaks = np.nanmax(values, axis=1)
    persistence = np.nanmean(values >= POWER_THRESHOLD_DB, axis=0)
    band_ranges = range_km[mask]
    persistent_ranges = band_ranges[persistence >= 0.5]
    return {
        "range_km": [lower_km, upper_km],
        "threshold_db": POWER_THRESHOLD_DB,
        "median_profile_peak_db": float(np.nanmedian(profile_peaks)),
        "maximum_db": float(np.nanmax(values)),
        "profiles_with_peak_above_threshold_fraction": float(
            np.nanmean(profile_peaks >= POWER_THRESHOLD_DB)
        ),
        "cells_above_threshold_fraction": float(
            np.nanmean(values >= POWER_THRESHOLD_DB)
        ),
        "maximum_gate_persistence_fraction": float(np.nanmax(persistence)),
        "persistent_gate_ranges_km": persistent_ranges.round(2).tolist(),
    }


def automatic_detections(
    metadata_reader: drf.DigitalMetadataReader,
    start_unix: float,
    end_unix: float,
) -> tuple[np.ndarray, np.ndarray]:
    phasecal = winds.load_phasecal()
    coordinates = winds.antenna_coordinates_ecef()
    position_differences = [
        coordinates[0, :] - coordinates[2, :],
        coordinates[0, :] - coordinates[3, :],
        coordinates[2, :] - coordinates[3, :],
    ]

    detections = []
    wind_rows = []
    prior_profile = None
    processing_start = start_unix - WIND_SEED_S

    for block_start in np.arange(processing_start, end_unix, WIND_BLOCK_S):
        block_end = min(block_start + WIND_BLOCK_S, end_unix)
        block_detections = winds.extract_detections_for_interval(
            metadata_reader,
            phasecal,
            position_differences,
            block_start,
            block_end,
            prior_wind_profile=prior_profile,
        )
        block_winds = (
            winds.estimate_wind_profile(block_detections, block_start)
            if len(block_detections)
            else winds.empty_wind_block(block_start)
        )
        wind_rows.append(block_winds)
        prior_profile = winds.latest_accepted_wind_profile(
            np.vstack(wind_rows)
        )
        if block_start >= start_unix and len(block_detections):
            detections.append(block_detections)

    selected = (
        np.vstack(detections)
        if detections
        else np.empty((0, 21), dtype=np.float64)
    )
    selected = winds.deduplicate_detections(selected)
    requested_winds = np.vstack(wind_rows)[
        np.vstack(wind_rows)[:, 0] >= start_unix
    ]
    return selected, requested_winds


def detection_statistics(detections: np.ndarray) -> dict:
    if len(detections) == 0:
        return {"count": 0}

    snr_db = 10.0 * np.log10(np.maximum(detections[:, 5], 1e-20))
    result = {
        "count": int(len(detections)),
        "unique_two_second_times": int(len(np.unique(detections[:, 0]))),
        "time_span_seconds": float(
            np.max(detections[:, 0]) - np.min(detections[:, 0])
        ),
        "height_km": {
            "minimum": float(np.min(detections[:, 3])),
            "median": float(np.median(detections[:, 3])),
            "maximum": float(np.max(detections[:, 3])),
        },
        "beamformed_snr_db": {
            "minimum": float(np.min(snr_db)),
            "median": float(np.median(snr_db)),
            "maximum": float(np.max(snr_db)),
        },
        "coherence": {
            "minimum": float(np.min(detections[:, 11])),
            "median": float(np.median(detections[:, 11])),
            "maximum": float(np.max(detections[:, 11])),
        },
        "radial_velocity_ms": {
            "minimum": float(np.min(detections[:, 4])),
            "median": float(np.median(detections[:, 4])),
            "maximum": float(np.max(detections[:, 4])),
        },
        "east_km": [
            float(np.min(detections[:, 1])),
            float(np.max(detections[:, 1])),
        ],
        "north_km": [
            float(np.min(detections[:, 2])),
            float(np.max(detections[:, 2])),
        ],
    }
    result["altitude_cut_counts"] = {
        f"{altitude:.0f}_plus_minus_2_km": int(
            np.count_nonzero(np.abs(detections[:, 3] - altitude) <= 2.0)
        )
        for altitude in winds.ALTITUDE_CUT_CENTERS_KM
    }
    return result


def wind_statistics(wind_rows: np.ndarray) -> dict:
    finite = np.isfinite(wind_rows[:, 2]) & np.isfinite(wind_rows[:, 3])
    return {
        "rows": int(len(wind_rows)),
        "finite_rows": int(np.count_nonzero(finite)),
        "finite_fraction": float(np.mean(finite)) if len(finite) else 0.0,
        "zonal_ms": (
            [
                float(np.min(wind_rows[finite, 2])),
                float(np.median(wind_rows[finite, 2])),
                float(np.max(wind_rows[finite, 2])),
            ]
            if np.any(finite)
            else None
        ),
        "meridional_ms": (
            [
                float(np.min(wind_rows[finite, 3])),
                float(np.median(wind_rows[finite, 3])),
                float(np.max(wind_rows[finite, 3])),
            ]
            if np.any(finite)
            else None
        ),
    }


def image_extent(
    times: np.ndarray,
    range_km: np.ndarray,
) -> tuple[float, float, float, float]:
    date_numbers = mdates.date2num(times.astype("datetime64[s]"))
    half_time_step = 0.5 * float(np.median(np.diff(date_numbers)))
    half_range_step = 0.5 * float(np.median(np.diff(range_km)))
    return (
        date_numbers[0] - half_time_step,
        date_numbers[-1] + half_time_step,
        range_km[0] - half_range_step,
        range_km[-1] + half_range_step,
    )


def plot_analysis(
    times: np.ndarray,
    range_km: np.ndarray,
    power_ratio_db: np.ndarray,
    detections: np.ndarray,
    start: dt.datetime,
    end: dt.datetime,
    output_path: Path,
) -> None:
    figure, axes = plt.subplots(
        3,
        1,
        figsize=(15, 13),
        facecolor="#070b14",
        constrained_layout=True,
    )
    for axis in axes:
        axis.set_facecolor("#050810")
        axis.tick_params(colors="#8fa1ba")
        axis.xaxis.label.set_color("#b7c5d9")
        axis.yaxis.label.set_color("#b7c5d9")
        axis.title.set_color("#edf4ff")
        for spine in axis.spines.values():
            spine.set_color("#314563")

    for axis, upper_range, limits, title in (
        (axes[0], 1500.0, (-3.0, 35.0), "Full-range coherent RTI"),
        (axes[1], 200.0, (-10.0, 20.0), "Mesospheric coherent RTI"),
    ):
        mask = (range_km >= 0) & (range_km <= upper_range)
        image = axis.imshow(
            power_ratio_db[:, mask].T,
            extent=image_extent(times, range_km[mask]),
            origin="lower",
            aspect="auto",
            interpolation="nearest",
            cmap="plasma",
            vmin=limits[0],
            vmax=limits[1],
        )
        axis.xaxis_date()
        axis.xaxis.set_major_formatter(
            mdates.DateFormatter("%H:%M", tz=dt.timezone.utc)
        )
        axis.set_ylim(0, upper_range)
        axis.set_ylabel("One-way range (km)")
        axis.set_title(title)
        colorbar = figure.colorbar(image, ax=axis, pad=0.01)
        colorbar.set_label("Power / 30–50 km background (dB)")

    if len(detections):
        axes[1].scatter(
            [
                dt.datetime.fromtimestamp(value, tz=dt.timezone.utc)
                for value in detections[:, 0]
            ],
            detections[:, 9],
            s=8,
            facecolors="none",
            edgecolors="#ffffff",
            linewidths=0.45,
            alpha=0.8,
        )

        velocity_limit = float(
            np.clip(np.percentile(np.abs(detections[:, 4]), 98), 20, 300)
        )
        positions = axes[2].scatter(
            detections[:, 1],
            detections[:, 2],
            c=detections[:, 4],
            s=10,
            cmap="seismic",
            vmin=-velocity_limit,
            vmax=velocity_limit,
            linewidths=0,
        )
        figure.colorbar(positions, ax=axes[2], pad=0.01).set_label(
            "Radial velocity (m/s)"
        )
    else:
        axes[2].text(
            0.5,
            0.5,
            "No automatic wind detections",
            transform=axes[2].transAxes,
            ha="center",
            va="center",
            color="#edf4ff",
        )

    axes[2].scatter([0], [0], marker="+", s=90, color="#edf4ff")
    axes[2].set_xlim(-200, 200)
    axes[2].set_ylim(-200, 200)
    axes[2].set_aspect("equal", adjustable="box")
    axes[2].set_xlabel("East–west position (km)")
    axes[2].set_ylabel("North–south position (km)")
    axes[2].set_title("Automatically accepted echo positions")
    axes[1].set_xlabel("Time (UTC)")
    axes[0].set_xlabel("Time (UTC)")
    figure.suptitle(
        f"Ramfjordmoen MF radar interval analysis\n"
        f"{start:%Y-%m-%d %H:%M}–{end:%H:%M} UTC",
        color="#edf4ff",
        fontsize=17,
    )
    temporary = output_path.with_suffix(".tmp.png")
    figure.savefig(temporary, dpi=160, facecolor=figure.get_facecolor())
    plt.close(figure)
    os.replace(temporary, output_path)


def main() -> None:
    args = parse_args()
    if args.end <= args.start:
        raise ValueError("end must be after start")

    start_unix = int(args.start.timestamp())
    end_unix = int(args.end.timestamp())
    times = np.arange(start_unix, end_unix, args.cadence, dtype=np.int64)
    if len(times) < 2:
        raise ValueError("interval must contain at least two cadence bins")

    raw_reader = drf.DigitalRFReader(args.raw_dir)
    power = np.asarray(
        [
            monitor.process_time_bin(
                raw_reader,
                int(timestamp),
                list(monitor.DEFAULT_CHANNELS),
                args.averages,
            )
            for timestamp in times
        ],
        dtype=np.float32,
    )
    range_km = monitor.range_vector_km(power.shape[1])
    background_mask = (
        (range_km >= monitor.THIRTY_MINUTE_NOISE_RANGE_KM[0])
        & (range_km <= monitor.THIRTY_MINUTE_NOISE_RANGE_KM[1])
    )
    background_power = float(np.mean(power[:, background_mask]))
    power_ratio_db = 10.0 * np.log10(
        np.maximum(power / background_power, 1e-20)
    )

    metadata_reader = drf.DigitalMetadataReader(args.metadata_dir)
    detections, wind_rows = automatic_detections(
        metadata_reader,
        start_unix,
        end_unix,
    )

    tag = f"{args.start:%Y%m%dT%H%M}_{args.end:%H%M}"
    output_dir = Path(args.output_root) / tag
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = output_dir / "interval_analysis.png"
    plot_analysis(
        times,
        range_km,
        power_ratio_db,
        detections,
        args.start,
        args.end,
        plot_path,
    )
    np.save(output_dir / "detections.npy", detections)
    np.save(output_dir / "winds.npy", wind_rows)

    metrics = {
        "start": args.start.isoformat().replace("+00:00", "Z"),
        "end": args.end.isoformat().replace("+00:00", "Z"),
        "background_method": "mean_power_over_complete_interval_and_30_to_50_km",
        "background_power": background_power,
        "mesospheric_slant_range": band_statistics(
            power_ratio_db,
            range_km,
            60.0,
            200.0,
        ),
        "ionospheric_slant_range": band_statistics(
            power_ratio_db,
            range_km,
            200.0,
            1500.0,
        ),
        "automatic_detections": detection_statistics(detections),
        "wind_estimates": wind_statistics(wind_rows),
    }
    (output_dir / "metrics.json").write_text(
        json.dumps(metrics, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(metrics, indent=2, allow_nan=False))
    print(f"Plot: {plot_path}")


if __name__ == "__main__":
    main()
