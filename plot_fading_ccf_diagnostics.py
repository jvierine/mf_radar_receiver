#!/usr/bin/env python3
"""Plot complete spaced-antenna CCF curves for the latest 15 minutes."""

from __future__ import annotations

import argparse
import datetime as dt
import os
from pathlib import Path

import digital_rf as drf
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

import fading_wind
import mf_conf as mc


DEFAULT_OUTPUT = Path(
    "/data2/plots/monitor/latest_fading_ccf_diagnostics_15m.png"
)
WINDOW_S = 15 * 60
CHANNEL_NAMES = ("ch1", "ch3", "ch4")
PAIR_COLORS = ("#35b8ff", "#ff9a3c", "#c57bff")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-dir", default=mc.xc_dir)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    reader = drf.DigitalMetadataReader(args.metadata_dir)
    end_unix = float(reader.get_bounds()[1]) / 1e6
    start_unix = end_unix - WINDOW_S
    sample_times, power = fading_wind.load_altitude_power(
        reader,
        start_unix,
        end_unix,
    )
    order = np.argsort(sample_times)
    sample_times = sample_times[order]
    power = power[order]
    sample_interval_s = float(np.median(np.diff(sample_times)))
    altitude_indices = np.flatnonzero(
        np.isclose(np.mod(fading_wind.ALTITUDE_KM, 10.0), 0.0)
    )
    figure, axes = plt.subplots(
        4,
        2,
        figsize=(17, 16),
        sharex=True,
        sharey=True,
        facecolor="#070b14",
        constrained_layout=True,
    )
    pair_labels = [
        f"{CHANNEL_NAMES[first]} → {CHANNEL_NAMES[second]}"
        for first, second in fading_wind.PAIR_INDICES
    ]
    for axis, altitude_index in zip(axes.ravel(), altitude_indices):
        altitude_km = fading_wind.ALTITUDE_KM[altitude_index]
        values = power[:, :, altitude_index]
        raw = fading_wind.correlation_diagnostics(
            values,
            sample_interval_s,
            conditioned=False,
        )
        conditioned = fading_wind.correlation_diagnostics(
            values,
            sample_interval_s,
            conditioned=True,
        )
        fit = fading_wind.fit_fading_window_diagnostics(
            values,
            sample_interval_s,
        )
        lag_seconds = conditioned["lag_seconds"]
        for pair_index, (label, color) in enumerate(
            zip(pair_labels, PAIR_COLORS)
        ):
            axis.plot(
                raw["lag_seconds"],
                raw["cross_correlations"][pair_index],
                color=color,
                linewidth=0.8,
                alpha=0.22,
            )
            axis.plot(
                lag_seconds,
                conditioned["cross_correlations"][pair_index],
                color=color,
                linewidth=1.6,
                label=label,
            )
            axis.plot(
                lag_seconds,
                fit["model_curves"][pair_index + 1],
                color=color,
                linewidth=1.0,
                linestyle="--",
            )
            axis.axvline(
                conditioned["delay_s"][pair_index],
                color=color,
                linewidth=0.7,
                alpha=0.6,
            )
        candidate_speed = np.hypot(
            fit["candidate_zonal_wind_ms"],
            fit["candidate_meridional_wind_ms"],
        )
        decision = "accepted" if fit["valid"] else "rejected"
        axis.set_title(
            f"{altitude_km:.0f} km · {decision} · "
            f"candidate |U|={candidate_speed:.1f} m/s"
        )
        axis.text(
            0.01,
            0.03,
            "min CCF={:.2f}  min prominence={:.2f}  "
            "max peak error={:.2f} s".format(
                np.min(conditioned["peak_correlation"]),
                np.min(conditioned["peak_prominence"]),
                np.max(fit["model_peak_error_s"]),
            ),
            transform=axis.transAxes,
            color="#b7c5d9",
            fontsize=9,
        )
        axis.axhline(0.0, color="#52657d", linewidth=0.7)
        axis.set_xlim(-fading_wind.MAX_LAG_S, fading_wind.MAX_LAG_S)
        axis.set_ylim(-0.35, 1.05)
        axis.set_facecolor("#050810")
        axis.grid(color="#1c2a3a", linewidth=0.5, alpha=0.7)
        axis.tick_params(colors="#8fa1ba")
        axis.title.set_color("#edf4ff")
        axis.set_ylabel("Normalized correlation", color="#b7c5d9")
        for spine in axis.spines.values():
            spine.set_color("#314563")
    for axis in axes[-1]:
        axis.set_xlabel("Lag (s)", color="#b7c5d9")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        ncol=3,
        frameon=False,
        labelcolor="#dbe7f8",
        bbox_to_anchor=(0.5, 0.975),
    )
    interval = (
        dt.datetime.fromtimestamp(start_unix, tz=dt.timezone.utc)
        .strftime("%Y-%m-%d %H:%M")
        + "–"
        + dt.datetime.fromtimestamp(end_unix, tz=dt.timezone.utc)
        .strftime("%H:%M UTC")
    )
    figure.suptitle(
        "Ramfjordmoen MF radar · complete three-baseline CCF curves · "
        f"{interval}\n"
        "faint: unfiltered power · solid: conditioned data · "
        "dashed: fitted FCA model",
        color="#edf4ff",
        fontsize=17,
        weight="semibold",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(".tmp.png")
    figure.savefig(temporary, dpi=150, facecolor=figure.get_facecolor())
    plt.close(figure)
    os.replace(temporary, args.output)
    print(f"CCF diagnostics: {args.output}")


if __name__ == "__main__":
    main()
