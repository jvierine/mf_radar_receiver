"""Three-dipole interferometric angle-of-arrival helpers."""

import numpy as np
from scipy.ndimage import maximum_filter

import jcoord


NGRID = 200
MAX_ZENITH_ANGLE_DEG = 45.0
PEAK_FILTER_SIZE = 5

u = np.linspace(-1, 1, num=NGRID)
v = np.linspace(-1, 1, num=NGRID)
uu, vv = np.meshgrid(u, v)
ww = np.sqrt(np.maximum(1 - uu**2 - vv**2, 0.0))

max_uhor = np.sin(np.radians(MAX_ZENITH_ANGLE_DEG))
inside_cone = np.hypot(uu, vv) < max_uhor
uu_search = uu[inside_cone].flatten()
vv_search = vv[inside_cone].flatten()
ww_search = ww[inside_cone].flatten()

# Candidate ENU directions expressed in ECEF for the ECEF antenna baselines.
ecef = jcoord.enu2ecef(
    69.58204,
    19.22283,
    0,
    uu_search,
    vv_search,
    ww_search,
)


def _match_grid(xc, posdiff, wavelength):
    """Return normalized phase-match amplitude over the AoA search grid."""
    klen = 2.0 * np.pi / wavelength
    measured = np.exp(1j * np.angle(xc))
    phasors = np.empty((len(xc), len(uu_search)), dtype=np.complex64)
    for index, baseline in enumerate(posdiff):
        predicted = klen * (
            baseline[0] * ecef[0, :]
            + baseline[1] * ecef[1, :]
            + baseline[2] * ecef[2, :]
        )
        phasors[index] = measured[index] * np.exp(1j * predicted)
    return np.abs(np.mean(phasors, axis=0)), phasors


def aoa_candidates(xc, posdiff, wavelength=100, max_candidates=12):
    """
    Return distinct interferometric grating-lobe candidates.

    Columns are east, north, up direction cosine, RMS phase residual, and
    normalized phase-match amplitude. Candidates are ordered by match.
    """
    match, phasors = _match_grid(xc, posdiff, wavelength)

    # Reconstruct the irregular cone samples on the regular u-v grid solely
    # for non-maximum suppression.
    match_grid = np.full(uu.shape, -np.inf, dtype=np.float32)
    match_grid[inside_cone] = match
    local_maximum = match_grid == maximum_filter(
        match_grid, size=PEAK_FILTER_SIZE, mode="constant", cval=-np.inf
    )
    peak_flat = np.flatnonzero(local_maximum[inside_cone])
    if len(peak_flat) == 0:
        peak_flat = np.array([int(np.argmax(match))])
    peak_flat = peak_flat[np.argsort(match[peak_flat])[::-1]][:max_candidates]

    candidates = []
    for candidate_index in peak_flat:
        rotated = phasors[:, candidate_index]
        common_phase = np.angle(np.mean(rotated))
        residual = np.angle(rotated * np.exp(-1j * common_phase))
        phase_residual = float(np.sqrt(np.mean(residual**2)))
        candidates.append(
            [
                uu_search[candidate_index],
                vv_search[candidate_index],
                ww_search[candidate_index],
                phase_residual,
                match[candidate_index],
            ]
        )
    return np.asarray(candidates)


def aoa(xc, posdiff, wavelength=100, return_quality=False):
    """Return the strongest three-dipole AoA candidate."""
    candidate = aoa_candidates(
        xc, posdiff, wavelength=wavelength, max_candidates=1
    )[0]
    result = tuple(candidate[:3])
    if return_quality:
        return (*result, candidate[3])
    return result
