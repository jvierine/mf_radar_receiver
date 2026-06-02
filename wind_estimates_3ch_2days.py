import digital_rf as drf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import h5py
import os
import datetime as dt

import jcoord
import mf_conf as mc
import image_help as ih




OUTDIR = "/data2/plots/realtime_winds/3ch_2day_plots"
os.makedirs(OUTDIR, exist_ok=True)

# First date to process (inclusive, UTC midnight).
# The script produces one 2-day plot per 2-day window starting here,
# rolling forward until no more complete data is available.
PROCESS_START = "2025-12-03"

WIND_DT       = 10 * 60   # 10-minute wind blocks (seconds)
PLOT_DURATION = 2          # days per plot
HEIGHT_RES    = 1.5        # km
HEIGHT_MIN    = 70
HEIGHT_MAX    = 122

READ_DT       = 60         # metadata chunk size, seconds

# Detection thresholds
NOISER0       = 50
NOISER1       = 110
NOISE_FMIN    = 5
SNR_THRESH    = 20
RANGE_MIN     = 50
RANGE_MAX     = 200
HEIGHT_DET_MIN = 70
HEIGHT_DET_MAX = 122
DOPPLER_MAX_HZ = 2.0

DOPPLER_SIGN  = -1
DOPPLER_TO_MS = mc.wavelength / 2.0

MIN_POINTS_PER_HEIGHT_BIN = 8

USE_BACKGROUND_REJECTION = True
BG_STRONG_FRAC_MAX       = 0.05




def unix_to_datetime(t):
    return dt.datetime.utcfromtimestamp(float(t))


def day_to_unix(date_str):
    """Return UTC midnight unix timestamp for a YYYY-MM-DD string."""
    return dt.datetime.strptime(date_str, "%Y-%m-%d").replace(
        tzinfo=dt.timezone.utc).timestamp()


def load_phasecal():
    if not os.path.exists("phasecal.h5"):
        raise FileNotFoundError("phasecal.h5 not found")
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


def weighted_lstsq(A, y, weights=None):
    if weights is None:
        return np.linalg.lstsq(A, y, rcond=None)[0]
    w = np.sqrt(np.asarray(weights))
    return np.linalg.lstsq(A * w[:, None], y * w, rcond=None)[0]


def fit_horizontal_wind(det, min_points=8):
    """
    Least-squares horizontal wind from radial velocities.

    det columns:
        0  time unix   4  vr m/s    8  los_w
        1  x km        5  snr       9  range km
        2  y km        6  los_u     10 doppler Hz
        3  height km   7  los_v
    """
    if det.shape[0] < min_points:
        return np.nan, np.nan, det.shape[0]
    A       = det[:, [6, 7]]
    y       = det[:, 4]
    weights = np.clip(det[:, 5], 1, None)
    try:
        U, V = weighted_lstsq(A, y, weights)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, det.shape[0]
    return U, V, det.shape[0]




def extract_detections_for_interval(dmt, phasecal, pos_diffs, start_unix, end_unix):
    """
    Extract meteor detections from start_unix to end_unix.

    Returns array with columns:
        0  time unix   4  vr m/s    8  los_w
        1  x km        5  snr       9  range km
        2  y km        6  los_u     10 doppler Hz
        3  height km   7  los_v
    """
    detections = []
    t0_us = int(start_unix * 1e6)
    t1_us = int(end_unix   * 1e6)

    for read_start in np.arange(t0_us, t1_us, int(READ_DT * 1e6)):
        read_end   = min(read_start + int(READ_DT * 1e6), t1_us)
        t_mid_unix = 0.5 * (read_start + read_end) / 1e6

        dd = dmt.read(read_start, read_end)

        for k in dd.keys():
            RDI1 = dd[k]["rdi1"] * np.exp(1j * phasecal[0])
            RDI2 = dd[k]["rdi2"] * np.exp(1j * phasecal[1])
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
                np.mean(np.abs(RDI1[0:fi0,  ri0:ri1]) ** 2.0)
                + np.mean(np.abs(RDI1[fi1:-1, ri0:ri1]) ** 2.0)
            )
            if noise_pwr <= 0 or not np.isfinite(noise_pwr):
                continue

            snr    = (np.abs(RDI1) ** 2.0 - noise_pwr) / noise_pwr
            rr, ff = np.meshgrid(rvec, fvec)

            if USE_BACKGROUND_REJECTION:
                strong_fraction = np.mean(snr > SNR_THRESH, axis=0)
                good_gate       = strong_fraction < BG_STRONG_FRAC_MAX
                good_cell       = good_gate[None, :]
            else:
                good_cell = np.ones_like(snr, dtype=bool)

            idx = np.where(
                (snr > SNR_THRESH)
                & (rr > RANGE_MIN) & (rr < RANGE_MAX)
                & (np.abs(ff) < DOPPLER_MAX_HZ)
                & good_cell
            )
            if len(idx[0]) == 0:
                continue

            xc13 = (RDI1[idx] * np.conj(RDI3[idx])).flatten()
            xc14 = (RDI1[idx] * np.conj(RDI4[idx])).flatten()
            xc34 = (RDI3[idx] * np.conj(RDI4[idx])).flatten()

            rrs  = rr[idx].flatten()
            ffs  = ff[idx].flatten()
            snrs = snr[idx].flatten()

            for mi in range(len(ffs)):
                los_u, los_v, los_w = ih.aoa(
                    [xc13[mi], xc14[mi], xc34[mi]],
                    pos_diffs,
                    wavelength=mc.wavelength,
                )

                x = rrs[mi] * los_u
                y = rrs[mi] * los_v
                z = rrs[mi] * los_w

                if z < HEIGHT_DET_MIN or z > HEIGHT_DET_MAX:
                    continue

                vr = DOPPLER_SIGN * DOPPLER_TO_MS * ffs[mi]

                detections.append([
                    t_mid_unix, x, y, z,
                    vr, snrs[mi],
                    los_u, los_v, los_w,
                    rrs[mi], ffs[mi],
                ])

    if len(detections) == 0:
        return np.empty((0, 11))
    return np.asarray(detections)



def estimate_wind_profile(det_block, block_start_unix):
    h_edges   = np.arange(HEIGHT_MIN, HEIGHT_MAX + HEIGHT_RES, HEIGHT_RES)
    h_centers = 0.5 * (h_edges[:-1] + h_edges[1:])
    block_center = block_start_unix + WIND_DT / 2
    rows = []
    for hi, h0 in enumerate(h_centers):
        m = (det_block[:, 3] >= h_edges[hi]) & (det_block[:, 3] < h_edges[hi + 1])
        U, V, N = fit_horizontal_wind(det_block[m], min_points=MIN_POINTS_PER_HEIGHT_BIN)
        rows.append([block_center, h0, U, V, N])
    return np.asarray(rows)



def plot_window(wind_rows, window_start_unix, window_end_unix, plot_file):
    if wind_rows is None or wind_rows.shape[0] == 0:
        print("  No wind data for this window, skipping plot.")
        return

    times   = np.unique(wind_rows[:, 0]); times.sort()
    heights = np.unique(wind_rows[:, 1]); heights.sort()

    U = np.full((len(heights), len(times)), np.nan)
    V = np.full_like(U, np.nan)

    for row in wind_rows:
        t, h, u, v, _ = row
        ti = np.where(times   == t)[0][0]
        hi = np.where(heights == h)[0][0]
        U[hi, ti] = u
        V[hi, ti] = v

    time_edges = np.concatenate([
        [times[0] - WIND_DT / 2],
        0.5 * (times[:-1] + times[1:]),
        [times[-1] + WIND_DT / 2],
    ])
    height_edges = np.concatenate([
        [heights[0] - HEIGHT_RES / 2],
        0.5 * (heights[:-1] + heights[1:]),
        [heights[-1] + HEIGHT_RES / 2],
    ])
    time_edges_dt = [unix_to_datetime(t) for t in time_edges]

    fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
    for ax, (Z, title) in zip(axes, [(U, "Eastward Wind"), (V, "Northward Wind")]):
        pc = ax.pcolormesh(time_edges_dt, height_edges, Z,
                           cmap="seismic", vmin=-100, vmax=100, shading="auto")
        ax.set_title(title)
        ax.set_ylabel("Height (km)")
        ax.set_ylim(HEIGHT_MIN, HEIGHT_MAX)
        ax.grid(True, alpha=0.2)
        fig.colorbar(pc, ax=ax).set_label("m/s")

    axes[-1].set_xlabel("Time UT")
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))

    t0_dt = unix_to_datetime(window_start_unix)
    t1_dt = unix_to_datetime(window_end_unix)
    fig.suptitle(
        f"MF meteor winds (3-ch)\n"
        f"{t0_dt:%Y-%m-%d %H:%M} – {t1_dt:%Y-%m-%d %H:%M} UT"
    )
    fig.tight_layout()
    fig.savefig(plot_file, dpi=200)
    plt.close(fig)
    print(f"  Saved plot: {plot_file}")



def main():
    phasecal = load_phasecal()

    coord = antenna_coordinates_ecef()
    pos_diffs = [
        coord[0, :] - coord[2, :],
        coord[0, :] - coord[3, :],
        coord[2, :] - coord[3, :],
    ]

    dmt = drf.DigitalMetadataReader(mc.xc_dir)

    bounds           = dmt.get_bounds()
    newest_available = bounds[1] / 1e6   # unix seconds

    window_start = day_to_unix(PROCESS_START)
    window_dur   = PLOT_DURATION * 86400  # 2 days in seconds

    print(f"Processing from {PROCESS_START} until data runs out.")
    print(f"Newest available data: {unix_to_datetime(newest_available)} UT\n")

    while True:
        window_end = window_start + window_dur

        # Only process fully covered windows
        if window_end > newest_available:
            print("Next window extends beyond available data — stopping.")
            break

        start_dt = unix_to_datetime(window_start)
        end_dt   = unix_to_datetime(window_end)
        print(f"=== Window {start_dt:%Y-%m-%d} – {end_dt:%Y-%m-%d} ===")

        date_tag  = start_dt.strftime("%Y-%m-%d")
        det_file  = os.path.join(OUTDIR, f"{date_tag}_detections_2day_3ch.npy")
        wind_file = os.path.join(OUTDIR, f"{date_tag}_winds_2day_3ch.npy")
        plot_file = os.path.join(OUTDIR, f"{date_tag}_mf_meteor_winds_2day_3ch.png")

        # Load saved progress for this window if it exists
        saved_winds = np.load(wind_file) if os.path.exists(wind_file) else None
        saved_dets  = np.load(det_file)  if os.path.exists(det_file)  else None

        if saved_winds is not None and saved_winds.size > 0:
            last_center = np.nanmax(saved_winds[:, 0])
            next_block  = last_center + WIND_DT / 2
        else:
            next_block  = np.floor(window_start / WIND_DT) * WIND_DT
            saved_winds = None
            saved_dets  = None

        all_winds = saved_winds
        all_dets  = saved_dets

        # Process every 10-minute block in this window
        while next_block < window_end:
            block_start = next_block
            block_end   = min(block_start + WIND_DT, window_end)

            print(f"  Block {unix_to_datetime(block_start):%Y-%m-%d %H:%M} – "
                  f"{unix_to_datetime(block_end):%H:%M} UT", end="", flush=True)

            det_block = extract_detections_for_interval(
                dmt, phasecal, pos_diffs, block_start, block_end
            )
            print(f"  →  {det_block.shape[0]} detections")

            if det_block.shape[0] > 0:
                wind_block = estimate_wind_profile(det_block, block_start)
            else:
                h_edges   = np.arange(HEIGHT_MIN, HEIGHT_MAX + HEIGHT_RES, HEIGHT_RES)
                h_centers = 0.5 * (h_edges[:-1] + h_edges[1:])
                block_center = block_start + WIND_DT / 2
                wind_block = np.asarray(
                    [[block_center, h, np.nan, np.nan, 0] for h in h_centers]
                )

            if all_dets is None or all_dets.size == 0:
                all_dets = det_block
            elif det_block.size > 0:
                all_dets = np.vstack([all_dets, det_block])

            if all_winds is None or all_winds.size == 0:
                all_winds = wind_block
            else:
                all_winds = np.vstack([all_winds, wind_block])

            # Save progress after every block
            np.save(det_file,  all_dets)
            np.save(wind_file, all_winds)

            next_block += WIND_DT

        # Plot the completed window
        if all_winds is not None and all_winds.size > 0:
            mask = (
                (all_winds[:, 0] >= window_start)
                & (all_winds[:, 0] <  window_end)
            )
            plot_window(all_winds[mask], window_start, window_end, plot_file)
        else:
            print("  No wind data for this window.")

        # Advance to the next 2-day window
        window_start = window_end


if __name__ == "__main__":
    main()
