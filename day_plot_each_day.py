import argparse
import os
import datetime

import digital_rf as drf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import mf_conf as mc


parser = argparse.ArgumentParser(
    description="Generate one day plot for each complete UTC day in the DigitalRF metadata directory."
)

parser.add_argument(
    "--outdir",
    type=str,
    default="images",
    help="Output directory for plots."
)

args = parser.parse_args()

# ------------------------
# Parameters
# ------------------------
ch = "ch1"                    # Channel, kept here in case you need it later
n_t = 24 * 60                 # Number of time bins, 1 per minute
n_rg = 1000                   # Number of range bins
rgs = np.arange(n_rg) * 1.5   # Range vector in km

os.makedirs(args.outdir, exist_ok=True)

# ------------------------
# DigitalRF metadata reader
# ------------------------
dmt = drf.DigitalMetadataReader(mc.xc_dir)

# Metadata bounds are in microseconds
bounds = dmt.get_bounds()
oldest_unix = bounds[0] / 1e6
newest_unix = bounds[1] / 1e6

print(f"Metadata bounds:")
print(f"  oldest_unix = {oldest_unix}")
print(f"  newest_unix = {newest_unix}")
print(f"  oldest UTC  = {datetime.datetime.utcfromtimestamp(oldest_unix)}")
print(f"  newest UTC  = {datetime.datetime.utcfromtimestamp(newest_unix)}")

# ------------------------
# Find complete UTC days in metadata bounds
# ------------------------
seconds_per_day = 24 * 60 * 60

# First possible complete day starts at the next UTC midnight after oldest_unix
first_midnight = (
    int(oldest_unix // seconds_per_day) * seconds_per_day
)

if oldest_unix > first_midnight:
    first_complete_day_t0 = first_midnight + seconds_per_day
else:
    first_complete_day_t0 = first_midnight

# Last possible complete day starts at the UTC midnight before newest_unix,
# but only if that whole day is included.
last_midnight = int(newest_unix // seconds_per_day) * seconds_per_day

complete_day_starts = []

t0_day = first_complete_day_t0
while t0_day + seconds_per_day <= newest_unix:
    complete_day_starts.append(t0_day)
    t0_day += seconds_per_day

if len(complete_day_starts) == 0:
    raise RuntimeError(
        "No complete UTC days found in metadata bounds. "
        "The directory may only contain partial-day data."
    )

print(f"Found {len(complete_day_starts)} complete UTC day(s).")

# ------------------------
# Plot settings
# ------------------------
plt.style.use("default")
plt.rcParams.update({
    "font.size": 12,
    "font.family": "serif",
    "figure.figsize": (10, 4),
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "axes.linewidth": 1.2,
})


# ------------------------
# Process each complete day
# ------------------------
for day_idx, t0 in enumerate(complete_day_starts, start=1):
    t1 = t0 + seconds_per_day

    date_str = datetime.datetime.utcfromtimestamp(t0).strftime("%Y-%m-%d")

    print("")
    print(f"Processing complete day {day_idx}/{len(complete_day_starts)}: {date_str}")
    print(f"  t0 = {t0}")
    print(f"  t1 = {t1}")

    # ------------------------
    # Initialize data array
    # ------------------------
    S = np.zeros([n_t, n_rg], dtype=np.float64)
    tvec = []

    # ------------------------
    # Loop over time
    # ------------------------
    for i in range(n_t):
        print(f"Processing {date_str}, minute {i+1}/{n_t}")

        start_us = int((t0 + i * 60) * 1e6)
        end_us = int(start_us + 60 * 1e6)

        dd = dmt.read(start_us, end_us)

        tvec.append(np.datetime64(int(t0 + i * 60), "s"))

        nk = 0
        for k in dd.keys():
            rdi = (
                np.abs(dd[k]["rdi1"])**2.0
                + np.abs(dd[k]["rdi2"])**2.0
                + np.abs(dd[k]["rdi3"])**2.0
                + np.abs(dd[k]["rdi4"])**2.0
            )

            S[i, :] += np.max(rdi, axis=0)
            nk += 1.0

        # Avoid divide by zero
        if nk > 0:
            S[i, :] /= nk
        else:
            S[i, :] = np.nan

    tvec = np.array(tvec)

    # ------------------------
    # Plot full range
    # ------------------------
    fig, ax = plt.subplots()

    dB = 10.0 * np.log10(S.T)

    pcm = ax.pcolormesh(
        tvec,
        rgs,
        dB,
        cmap="plasma",
        shading="auto",
        vmin=np.nanpercentile(dB, 10),
        vmax=np.nanpercentile(dB, 99),
    )

    ax.set_ylabel("Range (km)")
    ax.set_xlabel("Time (UTC)")
    ax.set_title(f"Range-Time Intensity for {date_str}")

    cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
    cbar.set_label("Received Power (dB)")

    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    fig.autofmt_xdate(rotation=0, ha="center")

    fig.tight_layout()

    full_range_filename = os.path.join(
        args.outdir,
        f"range_time_intensity_{date_str}.png"
    )

    fig.savefig(
        full_range_filename,
        dpi=600,
        bbox_inches="tight",
        transparent=False,
    )

    # ------------------------
    # Plot limited range
    # ------------------------
    ax.set_ylim([50, 200])
    fig.tight_layout()

    limited_range_filename = os.path.join(
        args.outdir,
        f"range_time_intensity_{date_str}_50_200km.png"
    )

    fig.savefig(
        limited_range_filename,
        dpi=600,
        bbox_inches="tight",
        transparent=False,
    )

    plt.close(fig)

    print(f"Saved:")
    print(f"  {full_range_filename}")
    print(f"  {limited_range_filename}")
