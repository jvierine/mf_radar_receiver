#!/usr/bin/env python3
"""Plot all candidates retained by the realtime Doppler–AoA search."""

from __future__ import annotations

import os
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np


VELOCITY_LIMIT_MS = 150.0
AMBIGUITY_PLOT = "latest_interferometry_ambiguities.png"


def column_map(handle):
    return {
        name: index
        for index, name in enumerate(
            handle.attrs["aoa_candidate_columns"].split(",")
        )
    }


def sorted_groups(candidates, columns):
    order = np.lexsort(
        (
            candidates[:, columns["range_index"]],
            candidates[:, columns["time_unix"]],
        )
    )
    candidates = candidates[order]
    keys = candidates[:, [columns["time_unix"], columns["range_index"]]]
    starts = np.r_[0, 1 + np.flatnonzero(np.any(np.diff(keys, axis=0), axis=1))]
    ends = np.r_[starts[1:], len(candidates)]
    return candidates, starts, ends


def style_axes(axes):
    for axis in np.ravel(axes):
        axis.set_facecolor("#050810")
        axis.tick_params(colors="#8fa1ba")
        axis.xaxis.label.set_color("#b7c5d9")
        axis.yaxis.label.set_color("#b7c5d9")
        axis.title.set_color("#edf4ff")
        for spine in axis.spines.values():
            spine.set_color("#314563")


def add_colorbar(figure, image, axis, label):
    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label(label, color="#b7c5d9")
    colorbar.ax.tick_params(colors="#8fa1ba")
    colorbar.outline.set_edgecolor("#314563")


def plot_ambiguities(candidates, columns, starts, ends, destination):
    counts = ends - starts
    representative = candidates[starts]
    sample_step = max(1, int(np.ceil(len(candidates) / 60_000)))
    sample = candidates[::sample_step]
    sample_range = sample[:, columns["range_km"]]
    east_km = sample_range * sample[:, columns["east_direction_cosine"]]
    north_km = sample_range * sample[:, columns["north_direction_cosine"]]

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(14, 10),
        facecolor="#070b14",
        constrained_layout=True,
    )
    style_axes(axes)
    time_values = [
        np.datetime64(int(value), "s").astype("datetime64[ms]").astype(object)
        for value in representative[:, columns["time_unix"]]
    ]
    count_plot = axes[0, 0].scatter(
        time_values,
        representative[:, columns["range_km"]],
        c=counts,
        s=5,
        cmap="viridis",
        vmin=1,
        vmax=12,
    )
    axes[0, 0].set_title("Retained joint Doppler–AoA candidates per cell")
    axes[0, 0].set_xlabel("Time (UTC)")
    axes[0, 0].set_ylabel("Round-trip range (km)")
    axes[0, 0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    add_colorbar(figure, count_plot, axes[0, 0], "Candidate count")

    position_plot = axes[0, 1].scatter(
        east_km,
        north_km,
        c=sample[:, columns["velocity_ms"]],
        s=2,
        alpha=0.35,
        cmap="seismic",
        vmin=-VELOCITY_LIMIT_MS,
        vmax=VELOCITY_LIMIT_MS,
        rasterized=True,
    )
    axes[0, 1].set_title("All retained ambiguous horizontal positions")
    axes[0, 1].set_xlabel("East (km)")
    axes[0, 1].set_ylabel("North (km)")
    axes[0, 1].set_aspect("equal", adjustable="box")
    add_colorbar(
        figure,
        position_plot,
        axes[0, 1],
        "Radial velocity (m/s)",
    )

    altitude_plot = axes[1, 0].scatter(
        sample_range,
        sample[:, columns["altitude_km"]],
        c=sample[:, columns["relative_power_db"]],
        s=2,
        alpha=0.4,
        cmap="magma",
        vmin=-6,
        vmax=0,
        rasterized=True,
    )
    axes[1, 0].set_title("WGS84 altitude allowed by each range")
    axes[1, 0].set_xlabel("Round-trip range (km)")
    axes[1, 0].set_ylabel("Altitude (km)")
    add_colorbar(figure, altitude_plot, axes[1, 0], "Relative power (dB)")

    axes[1, 1].hist(
        counts,
        bins=np.arange(0.5, 13.5, 1.0),
        color="#61d0ff",
        edgecolor="#070b14",
    )
    axes[1, 1].set_title("Ambiguity multiplicity")
    axes[1, 1].set_xlabel("Candidates retained per time–range cell")
    axes[1, 1].set_ylabel("Cells")
    figure.suptitle(
        "Three-dipole interferometry · ambiguity audit",
        color="#edf4ff",
        fontsize=19,
        weight="semibold",
    )
    temporary = destination.with_suffix(".tmp.png")
    figure.savefig(
        temporary,
        dpi=160,
        facecolor=figure.get_facecolor(),
    )
    plt.close(figure)
    os.replace(temporary, destination)


def main(data_file=None, plot_dir=None):
    data_file = Path(data_file or "/data2/plots/monitor/latest_doppler_15m.h5")
    plot_dir = Path(plot_dir or "/data2/plots/monitor")
    with h5py.File(data_file, "r") as handle:
        columns = column_map(handle)
        candidates = handle["aoa_candidates"][:]
    if not len(candidates):
        raise RuntimeError("realtime product contains no AoA candidates")
    candidates, starts, ends = sorted_groups(candidates, columns)
    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_ambiguities(
        candidates,
        columns,
        starts,
        ends,
        plot_dir / AMBIGUITY_PLOT,
    )
    print(
        f"Interferometry ambiguity audit: {len(candidates)} candidates"
        f" in {len(starts)} cells"
    )


if __name__ == "__main__":
    main()
