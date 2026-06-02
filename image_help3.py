
import numpy as n
import jcoord
import mf_conf as mc

ngrid = 300   # increased from 200 for better angular resolution
              # at 30 degree zenith limit with 1.52λ baselines

u = n.linspace(-1, 1, num=ngrid)
v = n.linspace(-1, 1, num=ngrid)
uu, vv = n.meshgrid(u, v)
ww = n.sqrt(n.maximum(1 - uu**2.0 - vv**2.0, 0))  # u²+v²+w²=1

max_zenith_angle = 25.0
max_uhor = n.sin(n.pi * max_zenith_angle / 180)
idx = n.where(n.sqrt(uu**2 + vv**2.0) < max_uhor)

uu_search = uu[idx].flatten()
vv_search = vv[idx].flatten()
ww_search = ww[idx].flatten()

# Convert direction cosines to ECEF unit vectors
ecef = jcoord.enu2ecef(69.58204, 19.22283, 0,
                        uu_search, vv_search, ww_search)


def aoa(xc, posdiff, wavelength=100):
    """
    Estimate angle of arrival from three cross-correlations.

    Parameters
    ----------
    xc       : list of 3 complex values [xc13, xc14, xc34]
    posdiff  : list of 3 ECEF baseline vectors (metres)
    wavelength: radar wavelength (metres)

    Returns
    -------
    (los_u, los_v, los_w) — direction cosines of best-fit direction
    """
    klen = 2.0 * n.pi / wavelength

    # Compute match function for each baseline independently
    # Shape: [n_search] for each
    match = n.ones(len(uu_search), dtype=n.float32)

    for i in range(len(xc)):
        # Phase of cross-correlation (the measured inter-antenna phase)
        phi_meas = n.angle(xc[i])

        # Predicted phase delay for each candidate direction
        phi_pred = klen * (
            posdiff[i][0] * ecef[0, :] +
            posdiff[i][1] * ecef[1, :] +
            posdiff[i][2] * ecef[2, :]
        )

        # Match function for this baseline: cos(phi_meas - phi_pred)
        # = 1 when prediction matches measurement, < 1 otherwise
        # Using cosine avoids the 2π ambiguity issue with phase differences
        m_i = n.cos(phi_meas - phi_pred)

        # Shift to [0, 1] range and multiply into combined match
        match *= (m_i + 1.0) / 2.0

    aoai = n.argmax(match)
    return uu_search[aoai], vv_search[aoai], ww_search[aoai]
