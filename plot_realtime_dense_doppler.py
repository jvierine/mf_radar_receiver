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

import mf_conf as mc
import plot_monitor_rti as monitor
from plot_dense_doppler import dense_doppler


WINDOW_S = 30 * 60
DISPLAY_LIMIT_MS = 200.0
PLOT_DIR = Path("/data2/plots/monitor")
COMBINED_PLOT = PLOT_DIR / "latest_snr_doppler_30m_0_200.png"
DATA_FILE = PLOT_DIR / "latest_doppler_30m.npz"
SNR_STATE_FILE = PLOT_DIR / "rti_30m_2s_state.npz"


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
    display_mask = (snr_ranges >= 0.0) & (snr_ranges <= 200.0)
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
        axis.set_ylim(0.0, 200.0)
        axis.set_ylabel("One-way range (km)")
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
        (0.0, 200.0),
    )
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
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
    print(f"Combined SNR/Doppler RTI: {COMBINED_PLOT}")
    print(f"Dense fit cells: {velocity.size}")


if __name__ == "__main__":
    main()
