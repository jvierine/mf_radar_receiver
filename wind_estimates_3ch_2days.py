import digital_rf as drf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import h5py
import os
import shutil
import datetime as dt
from scipy.ndimage import uniform_filter
from scipy.optimize import minimize_scalar

import jcoord
import mf_conf as mc
import image_help as ih




OUTDIR = "/data2/products/winds/3ch_coherent_2day"
os.makedirs(OUTDIR, exist_ok=True)
LATEST_SELECTED_RTI = "/data2/plots/monitor/latest_selected_wind_pixels.png"

# First date to process (inclusive, UTC midnight).
# The script produces one 2-day plot per 2-day window starting here,
# rolling forward until no more complete data is available.
PROCESS_START = "2026-07-22"

WIND_DT       = 10 * 60   # 10-minute wind blocks (seconds)
PLOT_DURATION = 2          # days per plot
HEIGHT_RES    = 1.5        # km
HEIGHT_MIN    = 60
HEIGHT_MAX    = 150

READ_DT       = 60         # metadata chunk size, seconds

# Detection thresholds
NOISER0       = 50
NOISER1       = 110
NOISE_FMIN    = 5
SNR_THRESH    = 20
RANGE_MIN     = 50
RANGE_MAX     = 200
HEIGHT_DET_MIN = 60
HEIGHT_DET_MAX = 150
DOPPLER_MAX_HZ = 2.0
COHERENCE_MIN  = 0.80
COHERENCE_WINDOW = (3, 3)  # Doppler bins x range gates
PHASE_CLOSURE_MAX = 0.35
AOA_PHASE_RESIDUAL_MAX = 0.35
WIND_CONTINUITY_SCALE_MS = 20.0
WIND_CONTINUITY_MAX_MS = 40.0

DOPPLER_SIGN  = -1
DOPPLER_TO_MS = mc.wavelength / 2.0
MAX_RADIAL_VELOCITY_MS = 300.0
MAX_FIT_DOPPLER_HZ = 2.0 * MAX_RADIAL_VELOCITY_MS / mc.wavelength
BEAMFORMED_SNR_MIN = SNR_THRESH

MIN_POINTS_PER_HEIGHT_BIN = 8
MIN_AZIMUTH_SECTORS = 3
AZIMUTH_SECTORS = 8
MAX_GEOMETRY_CONDITION = 10.0
MIN_VELOCITY_RESIDUAL_MS = 10.0

USE_BACKGROUND_REJECTION = True
BG_STRONG_FRAC_MAX       = 0.05

RADAR_LAT = 69.58204
RADAR_LON = 19.22283
RADAR_ALT_M = 0.0




def unix_to_datetime(t):
    return dt.datetime.fromtimestamp(float(t), tz=dt.timezone.utc)


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


def fit_complex_sinusoid(times, voltage):
    """
    Fit voltage = amplitude * exp(2j*pi*doppler*times) + residual.

    A dense coarse search prevents a bounded scalar optimizer from selecting
    the wrong local minimum. The hard radial-velocity bound also keeps the fit
    inside the unaliased region of the coherently integrated time series.
    """
    times = np.asarray(times, dtype=np.float64)
    voltage = np.asarray(voltage, dtype=np.complex128)
    frequencies = np.linspace(-MAX_FIT_DOPPLER_HZ, MAX_FIT_DOPPLER_HZ, 401)
    basis = np.exp(-2j * np.pi * frequencies[:, None] * times[None, :])
    coarse_score = np.abs(basis @ voltage) ** 2
    best = int(np.argmax(coarse_score))
    step = frequencies[1] - frequencies[0]
    lower = max(-MAX_FIT_DOPPLER_HZ, frequencies[best] - step)
    upper = min(MAX_FIT_DOPPLER_HZ, frequencies[best] + step)

    def residual_power(frequency):
        model = np.exp(2j * np.pi * frequency * times)
        amplitude = np.vdot(model, voltage) / np.vdot(model, model)
        residual = voltage - amplitude * model
        return float(np.mean(np.abs(residual) ** 2))

    result = minimize_scalar(
        residual_power,
        bounds=(lower, upper),
        method="bounded",
        options={"xatol": 1e-5},
    )
    frequency = float(result.x)
    model = np.exp(2j * np.pi * frequency * times)
    amplitude = np.vdot(model, voltage) / np.vdot(model, model)
    residual = voltage - amplitude * model
    noise_power = max(float(np.mean(np.abs(residual) ** 2)), 1e-30)
    snr = float(np.abs(amplitude) ** 2 / noise_power)
    return frequency, amplitude, snr


def beamform_three_dipoles(rti1, rti3, rti4, pos_diffs, east, north, up):
    """Steer ch1/ch3/ch4 to an AoA and add their complex voltages coherently."""
    direction_ecef = np.asarray(
        jcoord.enu2ecef(RADAR_LAT, RADAR_LON, RADAR_ALT_M, east, north, up)
    ).reshape(3)
    klen = 2.0 * np.pi / mc.wavelength
    phase_ch1 = klen * np.dot(pos_diffs[0], direction_ecef)
    phase_ch4 = klen * np.dot(-pos_diffs[2], direction_ecef)
    return (
        rti1 * np.exp(1j * phase_ch1)
        + rti3
        + rti4 * np.exp(1j * phase_ch4)
    ) / 3.0


def geographic_position(east_km, north_km, up_km):
    """Convert an ENU echo location relative to the radar to geodetic."""
    radar_ecef = np.asarray(
        jcoord.geodetic2ecef(RADAR_LAT, RADAR_LON, RADAR_ALT_M)
    ).reshape(3)
    offset_ecef = np.asarray(
        jcoord.enu2ecef(
            RADAR_LAT,
            RADAR_LON,
            RADAR_ALT_M,
            east_km * 1_000.0,
            north_km * 1_000.0,
            up_km * 1_000.0,
        )
    ).reshape(3)
    return jcoord.ecef2geodetic(*(radar_ecef + offset_ecef))


def weighted_lstsq(A, y, weights=None):
    if weights is None:
        return np.linalg.lstsq(A, y, rcond=None)[0]
    w = np.sqrt(np.asarray(weights))
    return np.linalg.lstsq(A * w[:, None], y * w, rcond=None)[0]


def local_mean(values, size=COHERENCE_WINDOW):
    """Locally average a real or complex Doppler-range array."""
    if np.iscomplexobj(values):
        return (
            uniform_filter(values.real, size=size, mode="nearest")
            + 1j * uniform_filter(values.imag, size=size, mode="nearest")
        )
    return uniform_filter(values, size=size, mode="nearest")


def three_dipole_coherence(rdi1, rdi3, rdi4):
    """
    Estimate normalized cross-spectral coherence on all dipole baselines.

    A single complex FFT pixel gives unit coherence by construction. The local
    3x3 average supplies nine Doppler-range looks and measures whether the
    inter-antenna phase is stable around the candidate echo.
    """
    p1 = local_mean(np.abs(rdi1) ** 2)
    p3 = local_mean(np.abs(rdi3) ** 2)
    p4 = local_mean(np.abs(rdi4) ** 2)

    xc13 = local_mean(rdi1 * np.conj(rdi3))
    xc14 = local_mean(rdi1 * np.conj(rdi4))
    xc34 = local_mean(rdi3 * np.conj(rdi4))

    tiny = np.finfo(np.float32).tiny
    coh13 = np.abs(xc13) / np.sqrt(np.maximum(p1 * p3, tiny))
    coh14 = np.abs(xc14) / np.sqrt(np.maximum(p1 * p4, tiny))
    coh34 = np.abs(xc34) / np.sqrt(np.maximum(p3 * p4, tiny))
    coherence = np.clip((coh13 + coh14 + coh34) / 3.0, 0.0, 1.0)

    return (xc13, xc14, xc34), coherence


def phase_closure(xc13, xc14, xc34):
    """Absolute ch1-ch3, ch1-ch4, ch3-ch4 phase-closure error."""
    return np.abs(np.angle(xc13 * np.conj(xc14) * xc34))


def prior_wind_at_height(prior_profile, height_km):
    """Return the nearest finite prior (zonal, meridional) wind estimate."""
    if prior_profile is None or len(prior_profile) == 0:
        return None
    finite = np.isfinite(prior_profile[:, 2]) & np.isfinite(prior_profile[:, 3])
    if not np.any(finite):
        return None
    rows = prior_profile[finite]
    nearest = rows[np.argmin(np.abs(rows[:, 1] - height_km))]
    return nearest[2], nearest[3]


def choose_continuous_aoa(
    xc,
    pos_diffs,
    range_km,
    radial_velocity_ms,
    prior_profile,
):
    """
    Resolve grating-lobe ambiguity using the preceding mean-wind profile.

    On startup, choose the strongest phase-consistent altitude-valid lobe.
    Afterwards, reject lobes whose measured radial velocity is inconsistent
    with the previous accepted zonal and meridional wind.
    """
    candidates = ih.aoa_candidates(
        xc, pos_diffs, wavelength=mc.wavelength
    )
    accepted = []
    for east, north, up, phase_residual, match in candidates:
        height_km = range_km * up
        if not HEIGHT_DET_MIN <= height_km <= HEIGHT_DET_MAX:
            continue
        if phase_residual > AOA_PHASE_RESIDUAL_MAX:
            continue

        prior = prior_wind_at_height(prior_profile, height_km)
        if prior is None:
            continuity_error = np.nan
            cost = (phase_residual / AOA_PHASE_RESIDUAL_MAX) ** 2 + (1.0 - match)
        else:
            predicted_velocity = prior[0] * east + prior[1] * north
            continuity_error = abs(radial_velocity_ms - predicted_velocity)
            if continuity_error > WIND_CONTINUITY_MAX_MS:
                continue
            cost = (
                (phase_residual / AOA_PHASE_RESIDUAL_MAX) ** 2
                + (continuity_error / WIND_CONTINUITY_SCALE_MS) ** 2
                + (1.0 - match)
            )
        accepted.append(
            (cost, east, north, up, phase_residual, match, continuity_error)
        )

    if not accepted:
        return None
    return min(accepted, key=lambda row: row[0])[1:]


def geometry_is_acceptable(A, weights):
    """Require azimuthal diversity and a well-conditioned wind inversion."""
    azimuth = np.mod(np.arctan2(A[:, 1], A[:, 0]), 2.0 * np.pi)
    sector_width = 2.0 * np.pi / AZIMUTH_SECTORS
    sectors = np.floor(azimuth / sector_width).astype(int)
    if len(np.unique(sectors)) < MIN_AZIMUTH_SECTORS:
        return False

    Aw = A * np.sqrt(weights)[:, None]
    singular_values = np.linalg.svd(Aw, compute_uv=False)
    if len(singular_values) < 2 or singular_values[-1] <= 0:
        return False
    return singular_values[0] / singular_values[-1] <= MAX_GEOMETRY_CONDITION


def fit_horizontal_wind(det, min_points=8):
    """
    Least-squares horizontal wind from radial velocities.

    det columns:
        0  time unix   4  vr m/s    8  los_w
        1  x km        5  snr       9  range km
        2  y km        6  los_u     10 doppler Hz
        3  height km   7  los_v       11 coherence
                                      12 phase closure
                                      13 AoA phase residual
                                      14 AoA match
                                      15 continuity residual
    """
    if det.shape[0] < min_points:
        return np.nan, np.nan, det.shape[0]
    A = det[:, [6, 7]]
    y = det[:, 4]
    weights = np.clip(det[:, 5], 1, None) * np.clip(det[:, 11], 0, 1) ** 2
    keep = np.ones(len(det), dtype=bool)

    for _ in range(3):
        count = np.count_nonzero(keep)
        if count < min_points or not geometry_is_acceptable(A[keep], weights[keep]):
            return np.nan, np.nan, count
        try:
            U, V = weighted_lstsq(A[keep], y[keep], weights[keep])
        except np.linalg.LinAlgError:
            return np.nan, np.nan, count

        residual = y - A @ np.array([U, V])
        center = np.median(residual[keep])
        mad = np.median(np.abs(residual[keep] - center))
        cutoff = max(MIN_VELOCITY_RESIDUAL_MS, 4.0 * 1.4826 * mad)
        updated = keep & (np.abs(residual - center) <= cutoff)
        if np.array_equal(updated, keep):
            return U, V, count
        keep = updated

    count = np.count_nonzero(keep)
    if count < min_points or not geometry_is_acceptable(A[keep], weights[keep]):
        return np.nan, np.nan, count
    U, V = weighted_lstsq(A[keep], y[keep], weights[keep])
    return U, V, count




def extract_detections_for_interval(
    dmt,
    phasecal,
    pos_diffs,
    start_unix,
    end_unix,
    prior_wind_profile=None,
):
    """
    Extract high-coherence three-dipole detections from start_unix to end_unix.

    Returns array with columns:
        0  time unix   4  vr m/s    8  los_w
        1  x km        5  snr       9  range km
        2  y km        6  los_u     10 doppler Hz
        3  height km   7  los_v       11 coherence
                                      12 phase closure
                                      13 AoA phase residual
                                      14 AoA match
                                      15 continuity residual
    """
    detections = []
    t0_us = int(start_unix * 1e6)
    t1_us = int(end_unix   * 1e6)

    for read_start in np.arange(t0_us, t1_us, int(READ_DT * 1e6)):
        read_end   = min(read_start + int(READ_DT * 1e6), t1_us)
        t_mid_unix = 0.5 * (read_start + read_end) / 1e6

        dd = dmt.read(read_start, read_end)

        for k in dd.keys():
            if not all(
                name in dd[k]
                for name in ("rti1", "rti3", "rti4", "tvec")
            ):
                continue
            RDI1 = dd[k]["rdi1"] * np.exp(1j * phasecal[0])
            RDI3 = dd[k]["rdi3"] * np.exp(1j * phasecal[2])
            RDI4 = dd[k]["rdi4"] * np.exp(1j * phasecal[3])
            RTI1 = dd[k]["rti1"] * np.exp(1j * phasecal[0])
            RTI3 = dd[k]["rti3"] * np.exp(1j * phasecal[2])
            RTI4 = dd[k]["rti4"] * np.exp(1j * phasecal[3])

            rvec = dd[k]["rvec"]
            fvec = dd[k]["fvec"]
            tvec = dd[k]["tvec"]

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

            cross_spectra, coherence = three_dipole_coherence(RDI1, RDI3, RDI4)
            closure = phase_closure(*cross_spectra)
            candidate_idx = np.where(
                (snr > SNR_THRESH)
                & (rr > RANGE_MIN) & (rr < RANGE_MAX)
                & (np.abs(ff) <= MAX_FIT_DOPPLER_HZ)
                & good_cell
                & (coherence >= COHERENCE_MIN)
                & (closure <= PHASE_CLOSURE_MAX)
            )
            if len(candidate_idx[0]) == 0:
                continue

            # Retain only the strongest locator pixel at each range gate. The
            # final Doppler comes from the sinusoid fit, not this FFT bin.
            selected = []
            candidate_quality = (
                snr[candidate_idx] * coherence[candidate_idx] ** 2
            )
            for range_index in np.unique(candidate_idx[1]):
                positions = np.flatnonzero(candidate_idx[1] == range_index)
                selected.append(positions[np.argmax(candidate_quality[positions])])
            selected = np.asarray(selected, dtype=int)
            doppler_indices = candidate_idx[0][selected]
            range_indices = candidate_idx[1][selected]
            idx = (doppler_indices, range_indices)

            xc13 = cross_spectra[0][idx].flatten()
            xc14 = cross_spectra[1][idx].flatten()
            xc34 = cross_spectra[2][idx].flatten()

            rrs  = rr[idx].flatten()
            ffs  = ff[idx].flatten()
            snrs = snr[idx].flatten()
            coherences = coherence[idx].flatten()
            closures = closure[idx].flatten()

            for mi in range(len(ffs)):
                locator_vr = DOPPLER_SIGN * DOPPLER_TO_MS * ffs[mi]
                solution = choose_continuous_aoa(
                    [xc13[mi], xc14[mi], xc34[mi]],
                    pos_diffs,
                    rrs[mi],
                    locator_vr,
                    prior_wind_profile,
                )
                if solution is None:
                    continue
                (
                    los_u,
                    los_v,
                    los_w,
                    phase_residual,
                    aoa_match,
                    continuity_residual,
                ) = solution

                x = rrs[mi] * los_u
                y = rrs[mi] * los_v
                z = rrs[mi] * los_w
                beamformed = beamform_three_dipoles(
                    RTI1[:, range_indices[mi]],
                    RTI3[:, range_indices[mi]],
                    RTI4[:, range_indices[mi]],
                    pos_diffs,
                    los_u,
                    los_v,
                    los_w,
                )
                fitted_doppler, _, beamformed_snr = fit_complex_sinusoid(
                    tvec, beamformed
                )
                if beamformed_snr < BEAMFORMED_SNR_MIN:
                    continue
                vr = DOPPLER_SIGN * DOPPLER_TO_MS * fitted_doppler

                prior = prior_wind_at_height(prior_wind_profile, z)
                if prior is not None:
                    predicted_velocity = prior[0] * los_u + prior[1] * los_v
                    continuity_residual = abs(vr - predicted_velocity)
                    if continuity_residual > WIND_CONTINUITY_MAX_MS:
                        continue

                latitude, longitude, geographic_altitude_m = geographic_position(
                    x, y, z
                )

                detections.append([
                    float(k) / 1e6, x, y, z,
                    vr, beamformed_snr,
                    los_u, los_v, los_w,
                    rrs[mi], fitted_doppler, coherences[mi],
                    closures[mi], phase_residual, aoa_match,
                    continuity_residual,
                    latitude, longitude, geographic_altitude_m / 1_000.0,
                    range_indices[mi], snrs[mi],
                ])

    if len(detections) == 0:
        return np.empty((0, 21))
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


def latest_accepted_wind_profile(wind_rows):
    """Return the newest finite wind estimate available at each height."""
    if wind_rows is None or wind_rows.size == 0:
        return None
    rows = []
    for height in np.unique(wind_rows[:, 1]):
        candidates = wind_rows[
            (wind_rows[:, 1] == height)
            & np.isfinite(wind_rows[:, 2])
            & np.isfinite(wind_rows[:, 3])
        ]
        if len(candidates):
            rows.append(candidates[np.argmax(candidates[:, 0])])
    return np.asarray(rows) if rows else None


def plot_selected_rti(detections, window_start_unix, window_end_unix, plot_file):
    """Plot only the automatically accepted pixels used by the wind fit."""
    if detections is None or detections.size == 0:
        return
    mask = (
        (detections[:, 0] >= window_start_unix)
        & (detections[:, 0] < window_end_unix)
    )
    selected = detections[mask]
    if len(selected) == 0:
        return

    snr_db = 10.0 * np.log10(np.maximum(selected[:, 5], 1e-20))
    figure, axis = plt.subplots(figsize=(14, 5), facecolor="#070b14")
    axis.set_facecolor("#050810")
    pixels = axis.scatter(
        [unix_to_datetime(value) for value in selected[:, 0]],
        selected[:, 9],
        c=snr_db,
        s=5,
        marker="s",
        linewidths=0,
        cmap="plasma",
        vmin=10,
        vmax=30,
        rasterized=True,
    )
    axis.set_title("Automatically selected wind-processing pixels")
    axis.set_xlabel("Time UT")
    axis.set_ylabel("One-way slant range (km)")
    axis.set_ylim(RANGE_MIN, RANGE_MAX)
    axis.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    axis.tick_params(colors="#8fa1ba")
    axis.xaxis.label.set_color("#b7c5d9")
    axis.yaxis.label.set_color("#b7c5d9")
    axis.title.set_color("#edf4ff")
    for spine in axis.spines.values():
        spine.set_color("#314563")
    axis.grid(color="#ffffff", alpha=0.07)
    colorbar = figure.colorbar(pixels, ax=axis)
    colorbar.set_label("Beamformed sinusoid SNR (dB)", color="#b7c5d9")
    colorbar.ax.tick_params(colors="#8fa1ba")
    colorbar.outline.set_edgecolor("#314563")
    figure.tight_layout()
    figure.savefig(plot_file, dpi=180, facecolor=figure.get_facecolor())
    plt.close(figure)
    os.makedirs(os.path.dirname(LATEST_SELECTED_RTI), exist_ok=True)
    temporary = LATEST_SELECTED_RTI + ".tmp"
    shutil.copy2(plot_file, temporary)
    os.replace(temporary, LATEST_SELECTED_RTI)



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
        f"MF coherent-echo winds (three dipoles)\n"
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

    newest_complete_block = np.floor(newest_available / WIND_DT) * WIND_DT

    while window_start < newest_complete_block:
        window_end = window_start + window_dur
        processing_end = min(window_end, newest_complete_block)

        start_dt = unix_to_datetime(window_start)
        end_dt   = unix_to_datetime(window_end)
        print(f"=== Window {start_dt:%Y-%m-%d} – {end_dt:%Y-%m-%d} ===")

        date_tag  = start_dt.strftime("%Y-%m-%d")
        det_file  = os.path.join(OUTDIR, f"{date_tag}_detections_2day_3ch.npy")
        wind_file = os.path.join(OUTDIR, f"{date_tag}_winds_2day_3ch.npy")
        plot_file = os.path.join(
            OUTDIR, f"{date_tag}_mf_coherent_echo_winds_2day_3ch.png"
        )
        selected_rti_file = os.path.join(
            OUTDIR, f"{date_tag}_selected_wind_pixels_rti_2day_3ch.png"
        )

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
        while next_block < processing_end:
            block_start = next_block
            block_end   = min(block_start + WIND_DT, window_end)

            print(f"  Block {unix_to_datetime(block_start):%Y-%m-%d %H:%M} – "
                  f"{unix_to_datetime(block_end):%H:%M} UT", end="", flush=True)

            prior_profile = latest_accepted_wind_profile(all_winds)
            det_block = extract_detections_for_interval(
                dmt,
                phasecal,
                pos_diffs,
                block_start,
                block_end,
                prior_wind_profile=prior_profile,
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

        # Plot the available part of this window.
        if all_winds is not None and all_winds.size > 0:
            mask = (
                (all_winds[:, 0] >= window_start)
                & (all_winds[:, 0] <  window_end)
            )
            plot_window(all_winds[mask], window_start, processing_end, plot_file)
            plot_selected_rti(
                all_dets,
                window_start,
                processing_end,
                selected_rti_file,
            )
        else:
            print("  No wind data for this window.")

        if processing_end < window_end:
            break

        # Advance to the next complete 2-day window.
        window_start = window_end


if __name__ == "__main__":
    main()
