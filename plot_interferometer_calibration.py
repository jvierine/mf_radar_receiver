#!/usr/bin/env python3
"""Plot time-block phase estimates from a vertical-echo calibration HDF5 file."""

import argparse
import datetime as dt

import h5py
import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np


CHANNEL_LABELS = (
    "ch1 MF3 dipole (reference)",
    "ch2 loop (excluded)",
    "ch3 MF1 dipole",
    "ch4 MF2 dipole",
)
CHANNEL_COLORS = ("#444444", "#d62728", "#1f77b4", "#2ca02c")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("calibration", help="calibration HDF5 input")
    parser.add_argument("output", help="PNG output")
    return parser.parse_args()


def unwrap_finite(values):
    result = np.full_like(values, np.nan, dtype=np.float64)
    finite = np.isfinite(values)
    result[finite] = np.unwrap(values[finite])
    return result


def main():
    args = parse_args()
    with h5py.File(args.calibration, "r") as handle:
        times = np.asarray(handle["block_start_unix"][()])
        phases = np.asarray(handle["block_phasecal_rad"][()])
        counts = np.asarray(handle["block_selected_cells"][()])
        phasecal = np.asarray(handle["phasecal"][()])
        start = float(handle.attrs["interval_start_unix"])
        end = float(handle.attrs["interval_end_unix"])
        doppler_max = float(handle.attrs["doppler_abs_max_hz"])

    time_values = [
        dt.datetime.fromtimestamp(value, tz=dt.timezone.utc) for value in times
    ]
    figure, (phase_axis, count_axis) = plt.subplots(
        2,
        1,
        figsize=(11, 6.8),
        sharex=True,
        gridspec_kw={"height_ratios": (3.2, 1.0)},
        constrained_layout=True,
    )

    for channel_index in range(1, 4):
        values = unwrap_finite(phases[:, channel_index])
        phase_axis.plot(
            time_values,
            values,
            marker="o",
            markersize=3.5,
            linewidth=1.4,
            color=CHANNEL_COLORS[channel_index],
            label=CHANNEL_LABELS[channel_index],
        )
        if channel_index in (2, 3):
            reference = phasecal[channel_index]
            finite = np.isfinite(values)
            if np.any(finite):
                reference += 2.0 * np.pi * np.round(
                    (np.nanmedian(values) - reference) / (2.0 * np.pi)
                )
            phase_axis.axhline(
                reference,
                color=CHANNEL_COLORS[channel_index],
                linewidth=1.0,
                linestyle="--",
                alpha=0.75,
            )

    phase_axis.set_ylabel("Phase correction relative to ch1 (rad)")
    phase_axis.set_title(
        "Vertical-echo phase calibration: dipoles stable, loop phase drifting"
    )
    phase_axis.grid(True, color="#d0d0d0", linewidth=0.6, alpha=0.7)
    phase_axis.legend(loc="best", frameon=False, ncol=1)

    count_axis.bar(
        time_values,
        counts,
        width=0.012,
        color="#777777",
        alpha=0.8,
    )
    count_axis.set_ylabel("Selected\nfits")
    count_axis.set_xlabel("Time (UTC)")
    count_axis.grid(True, axis="y", color="#d0d0d0", linewidth=0.6, alpha=0.7)
    count_axis.xaxis.set_major_locator(mdates.HourLocator(interval=1))
    count_axis.xaxis.set_major_formatter(
        mdates.DateFormatter("%H:%M", tz=dt.timezone.utc)
    )

    start_text = dt.datetime.fromtimestamp(
        start, tz=dt.timezone.utc
    ).strftime("%Y-%m-%d %H:%M")
    end_text = dt.datetime.fromtimestamp(
        end, tz=dt.timezone.utc
    ).strftime("%H:%M")
    figure.suptitle(
        f"{start_text}–{end_text} UTC · 10 s fits · |fD| ≤ {doppler_max:g} Hz",
        fontsize=11,
    )
    figure.savefig(args.output, dpi=180, facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    main()
