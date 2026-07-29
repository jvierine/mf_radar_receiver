
import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import digital_rf as drf
import logging
import os
import time

import mf_conf as mc
import calc_rti as crti


# ── Settings ──────────────────────────────────────────────────────────────────
BLOCK_S      = 2            # processing block length, seconds
N_SAMPLES    = BLOCK_S * 1_000_000   # raw samples per block at 1 MHz
CHANNELS     = ["ch1", "ch2", "ch3", "ch4"]
PLOT_DIR     = "/data2/plots"
SLEEP_S      = 0.5          # poll interval when waiting for new data, seconds

# ── Logging — errors and warnings only, written to file ───────────────────────
logging.basicConfig(
    filename="calc_xc.log",
    level=logging.WARNING,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def save_plots(rdis, rvec, tvec, i0):
    """Save per-channel RTI power plots for this block."""
    title = mc.unix2datestr(i0 / 1e6)
    dirname = mc.unix2dirname(i0 / 1e6)

    for ch, RDI in zip(CHANNELS, rdis):
        dB     = 10.0 * np.log10(np.abs(RDI.T) ** 2 + 1e-30)
        nfloor = np.nanmedian(dB)

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.pcolormesh(tvec, rvec, dB - nfloor,
                      cmap="plasma", vmin=-3, vmax=40, shading="auto")
        ax.set_title(f"{title}  [{ch}]")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Range (km)")
        ax.set_ylim(rvec[0], rvec[-1])
        cb = fig.colorbar(ax.collections[0], ax=ax)
        cb.set_label("Power (dB above noise)")
        fig.tight_layout()

        out_dir = os.path.join(PLOT_DIR, ch, dirname)
        os.makedirs(out_dir, exist_ok=True)
        fig.savefig(os.path.join(out_dir, f"rti-{int(i0/1e6):06d}.png"), dpi=100)
        plt.close(fig)


def process_block(d, dmw, i0):
    """
    Process one 2-second block across all channels.

    Returns True on success, False if the block could not be read.
    """
    rdis = []
    tvec = rvec = fvec = None

    for ch in CHANNELS:
        try:
            tvec, rvec, fvec, RTI, RDI = crti.rti(
                d, ch, i0,
                n_samples=N_SAMPLES,
                plot=False,
            )
        except Exception as e:
            logging.error("channel %s block %d: %s", ch, i0, e)
            return False
        rdis.append(RDI)

    # Write reduced data as Digital RF metadata
    data_out = {
        "rdi1": rdis[0],
        "rdi2": rdis[1],
        "rdi3": rdis[2],
        "rdi4": rdis[3],
        "rvec": rvec,
        "tvec": tvec,
        "fvec": fvec,
    }
    dmw.write(i0, data_out)

    # RTI plots use the RTI from the first channel read; re-read for plotting
    # using already-computed rdis (cheap: just the power image)
    save_plots(rdis, rvec, tvec, i0)

    return True


def main():
    os.makedirs(mc.xc_dir, exist_ok=True)

    d   = drf.DigitalRFReader(mc.raw_dir)
    dmw = drf.DigitalMetadataWriter(
        mc.xc_dir,
        3600,       # files per subdirectory (1 hour)
        BLOCK_S,    # samples per file (one block = BLOCK_S × 1e6 samples)
        1_000_000,  # sample rate numerator
        1,          # sample rate denominator
        "xc",
    )
    dmr = drf.DigitalMetadataReader(mc.xc_dir)

    # ── Determine where to start ──────────────────────────────────────────────
    b = d.get_bounds("ch1")
    try:
        db = dmr.get_bounds()
        i0 = db[1] + N_SAMPLES
    except Exception:
        i0 = b[0] + N_SAMPLES

    while True:
        b = d.get_bounds("ch1")

        # Skip blocks whose raw data has been overwritten by the ring buffer
        while i0 < b[0] + N_SAMPLES:
            logging.warning("block %d no longer in ring buffer, skipping", i0)
            i0 += N_SAMPLES

        # Wait until the full block is available
        while i0 + N_SAMPLES > b[1]:
            time.sleep(SLEEP_S)
            b = d.get_bounds("ch1")

        ok = process_block(d, dmw, i0)
        if not ok:
            logging.error("block %d failed, skipping", i0)

        i0 += N_SAMPLES


if __name__ == "__main__":
    main()
