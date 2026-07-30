
import numpy as np
import digital_rf as drf
import ctypes
import gc
import logging
import os
import sys
import time

import mf_conf as mc
import calc_rti as crti


# ── Settings ──────────────────────────────────────────────────────────────────
BLOCK_S      = 2            # processing block length, seconds
N_SAMPLES    = BLOCK_S * 1_000_000   # raw samples per block at 1 MHz
CHANNELS     = ["ch1", "ch2", "ch3", "ch4"]
SLEEP_S      = 0.5          # poll interval when waiting for new data, seconds
BACKFILL_LOAD_FRACTION = 0.75
BACKFILL_CURSOR_FILE = "/data2/metadata/xc_backfill_cursor.txt"
ALLOCATOR_TRIM_BLOCKS = 30
PROCESS_RECYCLE_S = 10 * 60

# ── Logging — errors and warnings only, written to file ───────────────────────
logging.basicConfig(
    filename="calc_xc.log",
    level=logging.WARNING,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def process_block(d, dmw, i0):
    """
    Process one 2-second block across all channels.

    Returns True on success, False if the block could not be read.
    """
    rtis = []
    rdis = []
    tvec = rvec = fvec = None

    for ch in CHANNELS:
        try:
            tvec, rvec, fvec, RTI, RDI = crti.rti(
                d, ch, i0,
                n_samples=N_SAMPLES,
            )
        except Exception as e:
            logging.error("channel %s block %d: %s", ch, i0, e)
            return False
        rtis.append(RTI)
        rdis.append(RDI)

    # Write reduced data as Digital RF metadata
    data_out = {
        "rdi1": rdis[0],
        "rdi2": rdis[1],
        "rdi3": rdis[2],
        "rdi4": rdis[3],
        "rti1": rtis[0],
        "rti2": rtis[1],
        "rti3": rtis[2],
        "rti4": rtis[3],
        "rvec": rvec,
        "tvec": tvec,
        "fvec": fvec,
    }
    dmw.write(i0, data_out)

    return True


def load_backfill_cursor(dmr, raw_start):
    try:
        with open(BACKFILL_CURSOR_FILE, encoding="ascii") as handle:
            return int(handle.read().strip())
    except (FileNotFoundError, ValueError):
        try:
            return int(dmr.get_bounds()[1] + N_SAMPLES)
        except Exception:
            return int(raw_start + N_SAMPLES)


def save_backfill_cursor(cursor):
    temporary = BACKFILL_CURSOR_FILE + ".tmp"
    with open(temporary, "w", encoding="ascii") as handle:
        handle.write(f"{cursor}\n")
    os.replace(temporary, BACKFILL_CURSOR_FILE)


def spare_cpu_available():
    cpu_count = os.cpu_count() or 1
    return os.getloadavg()[0] < cpu_count * BACKFILL_LOAD_FRACTION


def trim_allocator():
    """Return freed NumPy/HDF5 arenas to the OS when libc supports it."""
    gc.collect()
    try:
        libc = ctypes.CDLL(None)
        malloc_trim = libc.malloc_trim
        malloc_trim.argtypes = [ctypes.c_size_t]
        malloc_trim.restype = ctypes.c_int
        malloc_trim(0)
    except (AttributeError, OSError):
        pass


def recycle_process(reader):
    """Bound native-library retention by replacing this process in place."""
    try:
        reader.close()
    except Exception:
        pass
    trim_allocator()
    os.execv(sys.executable, [sys.executable, *sys.argv])


def latest_complete_block(bounds):
    """Start sample of the newest complete, offset-safe two-second block."""
    return (
        (bounds[1] - N_SAMPLES - crti.OFFSET) // N_SAMPLES
    ) * N_SAMPLES


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

    b = d.get_bounds("ch1")
    backfill_i0 = load_backfill_cursor(dmr, b[0])
    try:
        last_realtime_i0 = int(dmr.get_bounds()[1])
    except Exception:
        last_realtime_i0 = None
    process_started = time.monotonic()
    blocks_since_trim = 0

    while True:
        processed_block = False
        b = d.get_bounds("ch1")
        realtime_i0 = latest_complete_block(b)

        # The newest complete block always has priority.
        if realtime_i0 > b[0] and realtime_i0 != last_realtime_i0:
            if process_block(d, dmw, realtime_i0):
                last_realtime_i0 = realtime_i0
                processed_block = True
            else:
                logging.error("realtime block %d failed; retrying", realtime_i0)
                time.sleep(SLEEP_S)
                continue

        # Advance a persistent historical cursor only when raw data still
        # exists, it cannot collide with realtime, and the host has headroom.
        while backfill_i0 < b[0] + N_SAMPLES:
            logging.warning(
                "historical block %d was overwritten; advancing cursor",
                backfill_i0,
            )
            backfill_i0 += N_SAMPLES
            save_backfill_cursor(backfill_i0)

        if (
            backfill_i0 < realtime_i0 - N_SAMPLES
            and spare_cpu_available()
        ):
            if not process_block(d, dmw, backfill_i0):
                logging.error("historical block %d failed; skipping", backfill_i0)
            else:
                processed_block = True
            backfill_i0 += N_SAMPLES
            save_backfill_cursor(backfill_i0)
        else:
            time.sleep(SLEEP_S)

        if processed_block:
            blocks_since_trim += 1
            if blocks_since_trim >= ALLOCATOR_TRIM_BLOCKS:
                trim_allocator()
                blocks_since_trim = 0
        if time.monotonic() - process_started >= PROCESS_RECYCLE_S:
            recycle_process(d)


if __name__ == "__main__":
    main()
