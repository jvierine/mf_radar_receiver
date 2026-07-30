#!/usr/bin/env python3
"""Shared receiver geometry, calibration, and Doppler definitions."""

from __future__ import annotations

from pathlib import Path

import h5py
import jcoord
import numpy as np

import mf_conf as mc


MAX_RADIAL_VELOCITY_MS = 300.0
MAX_FIT_DOPPLER_HZ = 2.0 * MAX_RADIAL_VELOCITY_MS / mc.wavelength
DOPPLER_SIGN = -1.0
DOPPLER_TO_MS = mc.wavelength / 2.0
DOPPLER_FIT_DURATION_S = 1.0


def load_phasecal(path="phasecal.h5"):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")
    with h5py.File(path, "r") as handle:
        return np.copy(handle["phasecal"][()])


def antenna_coordinates_ecef():
    return np.asarray(
        [
            jcoord.geodetic2ecef(latitude, longitude, altitude)
            for latitude, longitude, altitude in mc.antenna_coords
        ],
        dtype=np.float64,
    )
