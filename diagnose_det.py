"""
diagnose_detections.py
======================
Reads one time interval from the xc metadata and plots the
distribution of every detection property WITHOUT any rejection.

This tells you where (if anywhere) layer returns and meteor echoes
are naturally separated in coherence, Doppler, zenith angle, SNR,
and height — so you can set thresholds based on actual data rather
than guessing.

Usage
-----
Edit START_UNIX / END_UNIX to a 10-minute window you know contains
both layer contamination and real meteors, then run:

    python3 diagnose_detections.py

Outputs
-------
  detection_scatter.png   — scatter plots: coherence vs Doppler,
                            coherence vs zenith angle, height vs Doppler
  detection_hist.png      — histograms of every property
  detection_data.h5       — raw detection array for further analysis
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import h5py
import digital_rf as drf
import mf_conf as mc
import image_help as ih
import jcoord
import datetime as dt
import os
from hdf5_store import save_array

# ============================================================
# Configure this block
# ============================================================

# Edit these to a window you want to inspect
# Leave END_UNIX = None to use START_UNIX + 10 minutes
import datetime
dt = datetime.datetime(2025, 12, 25, 14, 0, 0)
START_UNIX = dt.timestamp()    # e.g. 1724000000.0  — replace with your timestamp
END_UNIX   = None    # None -> START_UNIX + 600

# ============================================================
# Detection settings — match wind_estimates.py exactly
# ============================================================

NOISER0    = 50
NOISER1    = 110
NOISE_FMIN = 5
SNR_THRESH     = 20
RANGE_MIN      = 50
RANGE_MAX      = 120
HEIGHT_DET_MIN = 70
HEIGHT_DET_MAX = 120
DOPPLER_MAX_HZ = 2.0
READ_DT        = 60


# ============================================================
# Helpers
# ============================================================

def load_phasecal():
    with h5py.File("phasecal.h5", "r") as h:
        return np.copy(h["phasecal"][()])


def antenna_coordinates_ecef():
    coord = np.zeros((4, 3))
    for i in range(4):
        coord[i, :] = jcoord.geodetic2ecef(
            mc.antenna_coords[i][0],
            mc.antenna_coords[i][1],
            mc.antenna_coords[i][2],
        )
    return coord


# ============================================================
# Extract all detections with no rejection — record everything
# ============================================================

def extract_all(dmt, phasecal, pos_diffs, start_unix, end_unix):
    """
    Returns ndarray shape (N, 12):
        0  time_unix
        1  x km
        2  y km
        3  height km
        4  vr m/s
        5  snr
        6  los_u
        7  los_v
        8  los_w
        9  range km
        10 doppler Hz
        11 coherence (mean of 3 baselines)
    """
    rows = []
    t0_us = int(start_unix * 1e6)
    t1_us = int(end_unix   * 1e6)

    for read_start in np.arange(t0_us, t1_us, int(READ_DT * 1e6)):
        read_end   = min(read_start + int(READ_DT * 1e6), t1_us)
        t_mid_unix = 0.5 * (read_start + read_end) / 1e6

        dd = dmt.read(read_start, read_end)
        for k in dd.keys():
            RDI1 = dd[k]["rdi1"] * np.exp(1j * phasecal[0])
            RDI3 = dd[k]["rdi3"] * np.exp(1j * phasecal[2])
            RDI4 = dd[k]["rdi4"] * np.exp(1j * phasecal[3])
            rvec = dd[k]["rvec"]
            fvec = dd[k]["fvec"]

            try:
                ri0 = np.where(rvec > NOISER0)[0][0]
                ri1 = np.where(rvec > NOISER1)[0][0]
                fi0 = np.where(fvec < -NOISE_FMIN)[0][-1]
                fi1 = np.where(fvec >  NOISE_FMIN)[0][0]
            except IndexError:
                continue

            noise_pwr = 0.5 * (
                np.mean(np.abs(RDI1[0:fi0, ri0:ri1]) ** 2)
                + np.mean(np.abs(RDI1[fi1:-1, ri0:ri1]) ** 2)
            )
            if noise_pwr <= 0 or not np.isfinite(noise_pwr):
                continue

            snr    = (np.abs(RDI1) ** 2 - noise_pwr) / noise_pwr
            rr, ff = np.meshgrid(rvec, fvec)

            idx = np.where(
                (snr > SNR_THRESH)
                & (rr > RANGE_MIN) & (rr < RANGE_MAX)
                & (np.abs(ff) < DOPPLER_MAX_HZ)
            )
            if len(idx[0]) == 0:
                continue

            xc13 = (RDI1[idx] * np.conj(RDI3[idx])).flatten()
            xc14 = (RDI1[idx] * np.conj(RDI4[idx])).flatten()
            xc34 = (RDI3[idx] * np.conj(RDI4[idx])).flatten()

            p1  = np.abs(RDI1[idx]) ** 2 + 1e-30
            p3  = np.abs(RDI3[idx]) ** 2 + 1e-30
            p4  = np.abs(RDI4[idx]) ** 2 + 1e-30
            coh = (  np.abs(xc13) / np.sqrt(p1 * p3)
                   + np.abs(xc14) / np.sqrt(p1 * p4)
                   + np.abs(xc34) / np.sqrt(p3 * p4)) / 3.0

            rrs  = rr[idx].flatten()
            ffs  = ff[idx].flatten()
            snrs = snr[idx].flatten()

            for mi in range(len(ffs)):
                los_u, los_v, los_w = ih.aoa(
                    [xc13[mi], xc14[mi], xc34[mi]],
                    pos_diffs, wavelength=mc.wavelength,
                )
                z = rrs[mi] * los_w
                if z < HEIGHT_DET_MIN or z > HEIGHT_DET_MAX:
                    continue

                rows.append([
                    t_mid_unix,
                    rrs[mi] * los_u,
                    rrs[mi] * los_v,
                    z,
                    mc.wavelength / 2.0 * ffs[mi],
                    snrs[mi],
                    los_u, los_v, los_w,
                    rrs[mi],
                    ffs[mi],
                    coh[mi],
                ])

    return np.asarray(rows) if rows else np.empty((0, 12))


# ============================================================
# Plotting
# ============================================================

def plot_scatter(det):
    zenith_angle = np.degrees(np.arcsin(np.sqrt(det[:, 6]**2 + det[:, 7]**2)))
    coh          = det[:, 11]
    doppler      = det[:, 10]
    height       = det[:, 3]
    snr_db       = 10 * np.log10(np.clip(det[:, 5], 1e-6, None))

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle(f"Detection scatter  (N={len(det)})", fontsize=13)

    def sc(ax, x, y, xlabel, ylabel, title):
        h = ax.scatter(x, y, c=height, cmap="plasma",
                       s=4, alpha=0.5, rasterized=True,
                       norm=mcolors.Normalize(70, 120))
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        fig.colorbar(h, ax=ax, label="Height (km)")

    sc(axes[0, 0], coh,          doppler,
       "Coherence",               "Doppler (Hz)",
       "Coherence vs Doppler")

    sc(axes[0, 1], zenith_angle,  coh,
       "Zenith angle (deg)",      "Coherence",
       "Zenith angle vs Coherence")

    sc(axes[1, 0], height,        coh,
       "Height (km)",             "Coherence",
       "Height vs Coherence")

    sc(axes[1, 1], height,        doppler,
       "Height (km)",             "Doppler (Hz)",
       "Height vs Doppler")

    fig.tight_layout()
    fig.savefig("images/detection_scatter.png", dpi=150)
    plt.close(fig)
    print("Saved: detection_scatter.png")


def plot_histograms(det):
    zenith_angle = np.degrees(np.arcsin(np.sqrt(det[:, 6]**2 + det[:, 7]**2)))
    coh          = det[:, 11]
    doppler      = np.abs(det[:, 10])
    height       = det[:, 3]
    snr_db       = 10 * np.log10(np.clip(det[:, 5], 1e-6, None))

    props = [
        (coh,          "Coherence",          50, (0, 1)),
        (doppler,      "|Doppler| (Hz)",      50, (0, DOPPLER_MAX_HZ)),
        (height,       "Height (km)",         50, (HEIGHT_DET_MIN, HEIGHT_DET_MAX)),
        (zenith_angle, "Zenith angle (deg)",  50, (0, 90)),
        (snr_db,       "SNR (dB)",            50, None),
        (det[:, 9],    "Range (km)",          50, (RANGE_MIN, RANGE_MAX)),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle(f"Detection distributions  (N={len(det)})", fontsize=13)

    for ax, (data, label, bins, rng) in zip(axes.flat, props):
        ax.hist(data, bins=bins, range=rng, color="steelblue", edgecolor="none")
        ax.set_xlabel(label)
        ax.set_ylabel("Count")
        ax.set_title(label)
        # Mark percentiles
        for pct, col in [(25, "orange"), (50, "red"), (75, "orange")]:
            v = np.percentile(data, pct)
            ax.axvline(v, color=col, lw=1.5,
                       label=f"p{pct}={v:.3f}")
        ax.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig("images/detection_hist.png", dpi=150)
    plt.close(fig)
    print("Saved: detection_hist.png")


# ============================================================
# Main
# ============================================================

def main():
    phasecal  = load_phasecal()
    coord     = antenna_coordinates_ecef()
    pos_diffs = [coord[0]-coord[2], coord[0]-coord[3], coord[2]-coord[3]]

    dmt = drf.DigitalMetadataReader(mc.xc_dir)

    global START_UNIX, END_UNIX
    if START_UNIX is None:
        # Default: most recent 10-minute block in the metadata
        bounds     = dmt.get_bounds()
        END_UNIX   = bounds[1] / 1e6
        START_UNIX = END_UNIX - 600
        print(f"No START_UNIX set — using most recent 10 minutes: "
              f"{dt.datetime.utcfromtimestamp(START_UNIX)} – "
              f"{dt.datetime.utcfromtimestamp(END_UNIX)} UT")
    elif END_UNIX is None:
        END_UNIX = START_UNIX + 600

    print(f"Extracting detections (no rejection) ...")
    det = extract_all(dmt, phasecal, pos_diffs, START_UNIX, END_UNIX)
    print(f"Total detections: {len(det)}")

    if len(det) == 0:
        print("No detections found in this interval — try a different window.")
        return

    save_array("detection_data.h5", det)
    print("Saved: detection_data.h5")

    # Summary statistics
    print(f"\n--- Coherence ---")
    print(f"  min={det[:,11].min():.3f}  "
          f"p25={np.percentile(det[:,11],25):.3f}  "
          f"p50={np.percentile(det[:,11],50):.3f}  "
          f"p75={np.percentile(det[:,11],75):.3f}  "
          f"max={det[:,11].max():.3f}")

    print(f"--- |Doppler| (Hz) ---")
    print(f"  min={np.abs(det[:,10]).min():.3f}  "
          f"p50={np.percentile(np.abs(det[:,10]),50):.3f}  "
          f"max={np.abs(det[:,10]).max():.3f}")

    print(f"--- Height (km) ---")
    print(f"  min={det[:,3].min():.1f}  "
          f"p50={np.percentile(det[:,3],50):.1f}  "
          f"max={det[:,3].max():.1f}")

    zenith = np.degrees(np.arcsin(np.sqrt(det[:,6]**2 + det[:,7]**2)))
    print(f"--- Zenith angle (deg) ---")
    print(f"  min={zenith.min():.1f}  "
          f"p50={np.percentile(zenith,50):.1f}  "
          f"max={zenith.max():.1f}")

    plot_scatter(det)
    plot_histograms(det)
    print("\nDone. Look at detection_scatter.png and detection_hist.png")
    print("The coherence histogram will show you if there are two populations.")


if __name__ == "__main__":
    main()
