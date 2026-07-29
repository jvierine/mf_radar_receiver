#!/usr/bin/env python3
"""Generate dense fitted-Doppler RTIs for the latest 30 minutes."""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import digital_rf as drf
import numpy as np

import mf_conf as mc
from plot_dense_doppler import dense_doppler, plot_dense_map


WINDOW_S = 30 * 60
DISPLAY_LIMIT_MS = 200.0
PLOT_DIR = Path("/data2/plots/monitor")
FULL_PLOT = PLOT_DIR / "latest_doppler_30m_50_200.png"
FOCUSED_PLOT = PLOT_DIR / "latest_doppler_30m_75_125.png"
DATA_FILE = PLOT_DIR / "latest_doppler_30m.npz"


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
        (50.0, 200.0),
    )
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    plot_dense_map(
        times,
        ranges,
        velocity,
        start,
        end,
        FULL_PLOT,
        DISPLAY_LIMIT_MS,
    )
    focused = (ranges >= 75.0) & (ranges <= 125.0)
    plot_dense_map(
        times,
        ranges[focused],
        velocity[:, focused],
        start,
        end,
        FOCUSED_PLOT,
        DISPLAY_LIMIT_MS,
    )
    np.savez_compressed(
        DATA_FILE,
        time_unix=times,
        range_km=ranges,
        velocity_ms=velocity,
        sinusoid_snr=fit_snr,
    )
    print(f"Full Doppler RTI: {FULL_PLOT}")
    print(f"Focused Doppler RTI: {FOCUSED_PLOT}")
    print(f"Dense fit cells: {velocity.size}")


if __name__ == "__main__":
    main()
