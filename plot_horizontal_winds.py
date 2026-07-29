"""
plot_horizontal_winds.py
========================
For each 10-minute block and each 1.5 km height bin, produces a
horizontal wind vector scatter plot over a ±150 km field of view,
styled after Tsutsumi et al. (MST16, 2024) slide 20/21.

Directory structure:
    PLOT_BASE_DIR/
        YYYY-MM-DD/
            HH/
                YYYYMMDD_HHMM_h082.0km.png
                YYYYMMDD_HHMM_h083.5km.png
                ...

Usage:
    Called automatically from wind_estimate_4ch.py after each block,
    or run standalone to replot from saved detections_2days.npy.

Standalone:
    python3 plot_horizontal_winds.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import datetime as dt
import os

import mf_conf as mc


# ============================================================
# Settings — keep in sync with wind_estimate_4ch.py
# ============================================================

PLOT_BASE_DIR  = "horizontal_wind_plots"
DET_FILE       = "realtime_winds/4ch_snr_mask_detections_2days.npy"

HEIGHT_RES     = 10.0   # km  — aggregated for overview plots
HEIGHT_MIN     = 70
HEIGHT_MAX     = 122

WIND_DT        = 10 * 60   # seconds per block

FIELD_EXTENT   = 150        # km  — plot from -150 to +150 in x and y
DOPPLER_SIGN   = 1
DOPPLER_TO_MS  = mc.wavelength / 2.0

# Colormap limits for radial velocity
VMIN = -100   # m/s
VMAX =  100   # m/s


# ============================================================
# Directory helpers
# ============================================================

def _block_dir(block_start_unix):
    """
    Returns (date_str, hour_str, datetime_str) for a block start time.
    e.g. ("2025-12-03", "06", "20251203_0610")
    """
    d = dt.datetime.utcfromtimestamp(float(block_start_unix))
    date_str     = d.strftime("%Y-%m-%d")
    hour_str     = d.strftime("%H")
    datetime_str = d.strftime("%Y%m%d_%H%M")
    return date_str, hour_str, datetime_str


def _ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


# ============================================================
# Single height-bin plot
# ============================================================

def _plot_one_height(det_h, h0, block_start_unix, out_path):
    """
    Scatter plot of detections at one height bin.

    det_h : ndarray [N, 12]  — detections already filtered to this height bin
    h0    : float             — height bin centre (km)
    """
    fig, ax = plt.subplots(figsize=(6, 6))

    if det_h.shape[0] > 0:
        x  = det_h[:, 1]   # east-west km
        y  = det_h[:, 2]   # south-north km
        vr = det_h[:, 4]   # radial velocity m/s

        sc = ax.scatter(x, y, c=vr,
                        cmap="seismic",
                        vmin=VMIN, vmax=VMAX,
                        s=8, alpha=0.7, rasterized=True)
        cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label("Monostatic Doppler velocity (m/s)", fontsize=9)

    # Formatting
    ax.set_xlim(-FIELD_EXTENT, FIELD_EXTENT)
    ax.set_ylim(-FIELD_EXTENT, FIELD_EXTENT)
    ax.set_xlabel("West  ←  East (km)", fontsize=9)
    ax.set_ylabel("South  ←  North (km)", fontsize=9)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.25)
    ax.axhline(0, color="k", lw=0.5, alpha=0.4)
    ax.axvline(0, color="k", lw=0.5, alpha=0.4)

    t = dt.datetime.utcfromtimestamp(float(block_start_unix))
    ax.set_title(
        f"{t:%Y-%m-%d  %H:%M} – {(t + dt.timedelta(minutes=10)):%H:%M} UT\n"
        f"Height: {h0:.1f} km   N={det_h.shape[0]}",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


# ============================================================
# Main entry point — plot one full block
# ============================================================

def plot_block_horizontal(det_block, block_start_unix):
    """
    For every 1.5 km height bin, produce one horizontal scatter plot
    and save it into the appropriate directory.

    Parameters
    ----------
    det_block        : ndarray [N, 12]  — all detections for this block
    block_start_unix : float            — Unix timestamp of block start
    """
    h_edges   = np.arange(HEIGHT_MIN, HEIGHT_MAX + HEIGHT_RES, HEIGHT_RES)
    h_centers = 0.5 * (h_edges[:-1] + h_edges[1:])

    date_str, hour_str, datetime_str = _block_dir(block_start_unix)
    out_dir = _ensure_dir(
        os.path.join(PLOT_BASE_DIR, date_str, hour_str)
    )

    for hi, h0 in enumerate(h_centers):
        if det_block.shape[0] > 0:
            m     = (det_block[:, 3] >= h_edges[hi]) & \
                    (det_block[:, 3] <  h_edges[hi + 1])
            det_h = det_block[m]
        else:
            det_h = np.empty((0, 12))

        # Skip if no detections in this height bin
        if det_h.shape[0] == 0:
            continue

        fname    = f"{datetime_str}_h{h0:05.1f}km.png"
        out_path = os.path.join(out_dir, fname)
        _plot_one_height(det_h, h0, block_start_unix, out_path)

    print(f"  [horiz plots] saved {len(h_centers)} panels -> {out_dir}")


# ============================================================
# Standalone replot from saved detections
# ============================================================

if __name__ == "__main__":
    import sys

    if not os.path.exists(DET_FILE):
        print(f"No detections file found at {DET_FILE}")
        sys.exit(1)

    det_all = np.load(DET_FILE)
    print(f"Loaded {len(det_all)} detections from {DET_FILE}")

    # Group by 10-minute block
    block_times = np.floor(det_all[:, 0] / WIND_DT) * WIND_DT
    unique_blocks = np.unique(block_times)
    print(f"Found {len(unique_blocks)} blocks to plot")

    for block_start in unique_blocks:
        m         = block_times == block_start
        det_block = det_all[m]
        plot_block_horizontal(det_block, block_start)

    print("Done.")
