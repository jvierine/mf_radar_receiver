"""
image_help.py  —  Angle-of-Arrival estimation for the MF meteor radar.

Channel-to-antenna mapping:
    RDI1  →  ch1  (dipole,  E=143.91 m, N=-78.95 m  relative to ch3)
    RDI2  →  ch2  (loop,    E= 43.45 m, N= 42.25 m  relative to ch3)
    RDI3  →  ch3  (dipole,  reference, E=0, N=0)
    RDI4  →  ch4  (dipole,  E=  4.65 m, N=-164.57 m relative to ch3)

Strategy (hierarchical, two-stage):
  1. Disambiguation  — short baseline ch3-ch2 (60.6 m, d/λ = 0.56,
     unambiguous to ~63°) resolves which 2π-wrap is correct.
     Uses cross-correlation RDI3 × conj(RDI2).
  2. Fine estimate   — least-squares phase fit on the three dipole-only
     baselines (ch1-ch3, ch1-ch4, ch3-ch4), seeded from step 1.
     Uses RDI1×RDI3*, RDI1×RDI4*, RDI3×RDI4*.

The loop cross-correlation never enters the wind geometry.
"""

import numpy as np
from scipy.optimize import minimize

# ── Antenna ENU positions (metres, ch3 as reference) ────────────────────────
_ENU = {
    "ch1": np.array([ 143.91,  -78.95, 0.0]),   # dipole  (RDI1)
    "ch2": np.array([  43.45,   42.25, 0.0]),   # loop    (RDI2)
    "ch3": np.array([   0.00,    0.00, 0.0]),   # dipole  (RDI3, reference)
    "ch4": np.array([   4.65, -164.57, 0.0]),   # dipole  (RDI4)
}

# Dipole baselines [dE, dN]  —  b_ij = r_i - r_j
# These correspond to the cross-correlations computed in wind_estimates_test4.py:
#   xc13  =  RDI1 × conj(RDI3)   →  ch1 - ch3
#   xc14  =  RDI1 × conj(RDI4)   →  ch1 - ch4
#   xc34  =  RDI3 × conj(RDI4)   →  ch3 - ch4
_B_CH1_CH3 = _ENU["ch1"][:2] - _ENU["ch3"][:2]   # 164.1 m
_B_CH1_CH4 = _ENU["ch1"][:2] - _ENU["ch4"][:2]   # 163.5 m
_B_CH3_CH4 = _ENU["ch3"][:2] - _ENU["ch4"][:2]   # 164.6 m

_DIPOLE_BASELINES = np.array([_B_CH1_CH3, _B_CH1_CH4, _B_CH3_CH4])  # (3, 2)

# Short loop baseline for disambiguation
#   xc32  =  RDI3 × conj(RDI2)   →  ch3 - ch2  (60.6 m, d/λ = 0.56)
_B_CH3_CH2 = _ENU["ch3"][:2] - _ENU["ch2"][:2]   # 60.6 m

# Search cone
MAX_ZENITH_DEG = 45.0
_MAX_UHOR = np.sin(np.radians(MAX_ZENITH_DEG))


# ── Internal helpers ─────────────────────────────────────────────────────────

def predicted_phases(u, v, baselines, wavelength):
    """Phase predicted on each baseline for direction (u, v)."""
    klen = 2.0 * np.pi / wavelength
    return klen * (baselines[:, 0] * u + baselines[:, 1] * v)


def phase_cost(uv, measured_phases, baselines, wavelength):
    """Sum of squared wrapped residuals — continuous across 2π boundaries."""
    u, v = uv
    if u * u + v * v >= 1.0:
        return 1e9
    pred = predicted_phases(u, v, baselines, wavelength)
    residuals = np.angle(np.exp(1j * (measured_phases - pred)))
    return float(np.dot(residuals, residuals))


def initial_guess(measured_phases, baselines, wavelength):
    """Analytic 2x2 solve for a fast starting point."""
    klen = 2.0 * np.pi / wavelength
    try:
        uv0 = np.linalg.solve(baselines[:2], measured_phases[:2] / klen)
    except np.linalg.LinAlgError:
        return np.array([0.0, 0.0])
    norm = np.hypot(uv0[0], uv0[1])
    if norm >= 1.0:
        uv0 = uv0 / norm * 0.95
    return uv0


def optimise(measured_phases, baselines, wavelength, uv0):
    """Nelder-Mead minimisation of wrapped phase residuals."""
    res = minimize(
        phase_cost,
        uv0,
        args=(measured_phases, baselines, wavelength),
        method="Nelder-Mead",
        options={"xatol": 1e-7, "fatol": 1e-12, "maxiter": 600},
    )
    return res.x


def clamp_to_cone(uv, max_uhor):
    r = np.hypot(uv[0], uv[1])
    if r > max_uhor:
        return uv / r * max_uhor * 0.999
    return uv


def aoa(xc_dipole, xc_loop, wavelength):
    """
    Estimate angle of arrival.

    Parameters
    ----------
    xc_dipole : sequence of 3 complex
        Cross-correlations between dipole pairs at one (Doppler, range) pixel:
            xc_dipole[0]  =  RDI1 × conj(RDI3)   →  ch1-ch3 baseline
            xc_dipole[1]  =  RDI1 × conj(RDI4)   →  ch1-ch4 baseline
            xc_dipole[2]  =  RDI3 × conj(RDI4)   →  ch3-ch4 baseline

    xc_loop : complex
        RDI3 × conj(RDI2)  →  ch3-ch2 baseline (60.6 m).
        Used only for grating-lobe disambiguation; not used in the wind fit.

    wavelength : float
        Radar wavelength in metres.

    Returns
    -------
    u, v, w : float
        East, North, Up direction cosines.
    phase_residual : float
        RMS wrapped phase residual on the dipole baselines (radians).
        Use this to reject detections with inconsistent phase across baselines.
    """
    meas_dipole = np.array([np.angle(xc) for xc in xc_dipole])  # (3,)

    # Stage 1 — coarse direction from short loop baseline + one long baseline.
    # The loop baseline (d/λ = 0.56) is unambiguous to 63°, so this coarse
    # solution always lands in the correct grating-lobe cell.
    coarse_baselines = np.array([_B_CH3_CH2, _B_CH1_CH3])
    meas_coarse      = np.array([np.angle(xc_loop), meas_dipole[0]])

    uv0_coarse = initial_guess(meas_coarse, coarse_baselines, wavelength)
    uv_coarse  = optimise(meas_coarse, coarse_baselines, wavelength, uv0_coarse)
    uv_coarse  = clamp_to_cone(uv_coarse, _MAX_UHOR)

    # Stage 2 — fine direction from dipole baselines only, seeded by stage 1.
    uv_fine = optimise(meas_dipole, _DIPOLE_BASELINES, wavelength, uv_coarse)
    uv_fine = clamp_to_cone(uv_fine, _MAX_UHOR)

    u, v = uv_fine
    w = np.sqrt(max(0.0, 1.0 - u * u - v * v))

    # Phase residual — quality metric for the dipole-only fit
    pred      = predicted_phases(u, v, _DIPOLE_BASELINES, wavelength)
    residuals = np.angle(np.exp(1j * (meas_dipole - pred)))
    phase_residual = float(np.sqrt(np.mean(residuals ** 2)))

    return u, v, w, phase_residual
