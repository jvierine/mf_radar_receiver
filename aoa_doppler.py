#!/usr/bin/env python3
"""Joint ten-second Doppler and three-dipole interferometric AoA search."""

from __future__ import annotations

from dataclasses import dataclass

import jcoord
import numpy as np
from scipy.ndimage import maximum_filter, maximum_filter1d

import mf_conf as mc


RADAR_LAT = 69.58204
RADAR_LON = 19.22283
RADAR_ALT_M = 0.0
DIPOLE_CHANNEL_INDICES = np.asarray((0, 2, 3))
DIPOLE_NAMES = ("ch1", "ch3", "ch4")
REFERENCE_DIPOLE = 1  # ch3 within the complete four-channel phase vector

ALTITUDE_MIN_KM = 50.0
ALTITUDE_MAX_KM = 150.0
MAX_ZENITH_DEG = 70.0
GRID_SIZE = 101
DIRECTION_PEAK_FILTER_SIZE = 5
DOPPLER_PEAK_FILTER_SIZE = 5
MAX_AMBIGUITIES_PER_RANGE = 12
AMBIGUITY_RELATIVE_POWER_DB = -6.0

CANDIDATE_COLUMNS = (
    "range_index",
    "range_km",
    "frequency_hz",
    "velocity_ms",
    "east_direction_cosine",
    "north_direction_cosine",
    "up_direction_cosine",
    "latitude_deg",
    "longitude_deg",
    "altitude_km",
    "beam_power_ratio",
    "relative_power_db",
    "aoa_match",
)


def _ecef_altitude_m(x, y, z):
    """Vectorized form of jcoord.ecef2geodetic's WGS84 altitude."""
    radius = np.sqrt(x * x + y * y)
    ellipsoid_difference = jcoord.a * jcoord.a - jcoord.b * jcoord.b
    f_term = 54.0 * jcoord.b * jcoord.b * z * z
    g_term = (
        radius * radius
        + (1.0 - jcoord.esq) * z * z
        - jcoord.esq * ellipsoid_difference
    )
    c_term = (
        jcoord.esq
        * jcoord.esq
        * f_term
        * radius
        * radius
        / np.maximum(g_term**3, np.finfo(float).tiny)
    )
    s_term = np.cbrt(1.0 + c_term + np.sqrt(c_term * c_term + 2.0 * c_term))
    p_term = f_term / (
        3.0
        * (s_term + 1.0 / s_term + 1.0) ** 2
        * g_term
        * g_term
    )
    q_term = np.sqrt(1.0 + 2.0 * jcoord.esq * jcoord.esq * p_term)
    r0 = (
        -(p_term * jcoord.esq * radius) / (1.0 + q_term)
        + np.sqrt(
            0.5 * jcoord.a * jcoord.a * (1.0 + 1.0 / q_term)
            - p_term
            * (1.0 - jcoord.esq)
            * z
            * z
            / (q_term * (1.0 + q_term))
            - 0.5 * p_term * radius * radius
        )
    )
    u_term = np.sqrt((radius - jcoord.esq * r0) ** 2 + z * z)
    v_term = np.sqrt(
        (radius - jcoord.esq * r0) ** 2 + (1.0 - jcoord.esq) * z * z
    )
    return u_term * (1.0 - jcoord.b * jcoord.b / (jcoord.a * v_term))


@dataclass(frozen=True)
class DirectionGrid:
    axis: np.ndarray
    u: np.ndarray
    v: np.ndarray
    w: np.ndarray
    regular_flat_index: np.ndarray
    direction_ecef: np.ndarray
    steering: np.ndarray
    radar_ecef: np.ndarray

    @property
    def shape(self):
        return (len(self.axis), len(self.axis))

    def geometry_for_ranges(self, range_km):
        """Return feasible mask and WGS84 altitude for every range/direction."""
        range_km = np.asarray(range_km, dtype=np.float64)
        altitude = np.empty((len(range_km), len(self.u)), dtype=np.float32)
        for range_index, value in enumerate(range_km):
            target = (
                self.radar_ecef[:, None]
                + value * 1000.0 * self.direction_ecef
            )
            altitude[range_index] = (
                _ecef_altitude_m(target[0], target[1], target[2]) / 1000.0
            )
        feasible = (
            (altitude >= ALTITUDE_MIN_KM)
            & (altitude <= ALTITUDE_MAX_KM)
        )
        return feasible, altitude

    def geographic_position(self, range_km, direction_index):
        target = (
            self.radar_ecef
            + range_km * 1000.0 * self.direction_ecef[:, direction_index]
        )
        return np.asarray(jcoord.ecef2geodetic(*target), dtype=np.float64)


def antenna_coordinates_ecef():
    return np.asarray(
        [
            jcoord.geodetic2ecef(latitude, longitude, altitude)
            for latitude, longitude, altitude in mc.antenna_coords
        ],
        dtype=np.float64,
    )


def build_direction_grid(grid_size=GRID_SIZE):
    """Build the regular ENU direction-cosine grid and steering matrix."""
    if grid_size < 3 or grid_size % 2 == 0:
        raise ValueError("grid_size must be an odd integer of at least three")
    axis = np.linspace(-1.0, 1.0, grid_size)
    uu, vv = np.meshgrid(axis, axis)
    horizontal_limit = np.sin(np.radians(MAX_ZENITH_DEG))
    inside = np.hypot(uu, vv) <= horizontal_limit
    regular_flat_index = np.flatnonzero(inside)
    u = uu.flat[regular_flat_index]
    v = vv.flat[regular_flat_index]
    w = np.sqrt(np.maximum(0.0, 1.0 - u * u - v * v))
    direction_ecef = np.asarray(
        jcoord.enu2ecef(
            RADAR_LAT,
            RADAR_LON,
            RADAR_ALT_M,
            u,
            v,
            w,
        ),
        dtype=np.float64,
    ).reshape(3, -1)

    coordinates = antenna_coordinates_ecef()
    dipole_coordinates = coordinates[DIPOLE_CHANNEL_INDICES]
    reference = dipole_coordinates[1]
    baselines = dipole_coordinates - reference[None, :]
    phase = (
        2.0
        * np.pi
        / mc.wavelength
        * (baselines @ direction_ecef)
    )
    steering = np.exp(1j * phase.T)
    radar_ecef = np.asarray(
        jcoord.geodetic2ecef(RADAR_LAT, RADAR_LON, RADAR_ALT_M),
        dtype=np.float64,
    ).reshape(3)
    return DirectionGrid(
        axis=axis,
        u=u,
        v=v,
        w=w,
        regular_flat_index=regular_flat_index,
        direction_ecef=direction_ecef,
        steering=steering,
        radar_ecef=radar_ecef,
    )


def _top_doppler_indices(incoherent_power):
    """Return every local Doppler maximum, strongest upper bound first."""
    local = incoherent_power == maximum_filter1d(
        incoherent_power,
        size=DOPPLER_PEAK_FILTER_SIZE,
        mode="nearest",
    )
    indices = np.flatnonzero(local & np.isfinite(incoherent_power))
    if len(indices) == 0:
        return np.asarray([int(np.nanargmax(incoherent_power))])
    order = np.argsort(incoherent_power[indices])[::-1]
    return indices[order]


def _direction_local_maxima(grid, values, feasible):
    regular = np.full(grid.shape, -np.inf, dtype=np.float64)
    flat = regular.ravel()
    flat[grid.regular_flat_index[feasible]] = values
    local = regular == maximum_filter(
        regular,
        size=DIRECTION_PEAK_FILTER_SIZE,
        mode="constant",
        cval=-np.inf,
    )
    regular_peaks = np.flatnonzero(local.ravel() & np.isfinite(regular.ravel()))
    if len(regular_peaks) == 0:
        return np.asarray([int(np.nanargmax(values))])
    lookup = np.full(regular.size, -1, dtype=np.int64)
    lookup[grid.regular_flat_index] = np.arange(len(grid.u))
    return lookup[regular_peaks]


def _parabolic_frequency(frequency, power, center_index):
    if center_index <= 0 or center_index >= len(frequency) - 1:
        return float(frequency[center_index])
    lower, center, upper = power
    denominator = lower - 2.0 * center + upper
    if abs(denominator) <= 1e-30:
        return float(frequency[center_index])
    offset = np.clip(0.5 * (lower - upper) / denominator, -1.0, 1.0)
    return float(
        frequency[center_index]
        + offset * (frequency[center_index + 1] - frequency[center_index])
    )


def joint_aoa_doppler_search(
    times,
    dipole_voltage,
    range_km,
    noise_power,
    grid,
    maximum_frequency_hz,
):
    """
    Search Doppler and all feasible AoA maxima for one ten-second voltage group.

    Returns best-frequency, best-velocity, best beam-power ratio, and a numeric
    candidate table with columns defined by ``CANDIDATE_COLUMNS``.
    """
    times = np.asarray(times, dtype=np.float64)
    times = times - times[0]
    voltage = np.asarray(dipole_voltage, dtype=np.complex128)
    range_km = np.asarray(range_km, dtype=np.float64)
    noise_power = np.asarray(noise_power, dtype=np.float64)
    if voltage.shape != (len(times), 3, len(range_km)):
        raise ValueError("dipole_voltage must have shape (time, 3, range)")
    if noise_power.shape != (3,) or np.any(noise_power <= 0):
        raise ValueError("noise_power must contain three positive values")

    sample_interval = float(np.median(np.diff(times)))
    nfft = max(2048, 2 ** int(np.ceil(np.log2(len(times) * 8))))
    frequency_full = np.fft.fftshift(
        np.fft.fftfreq(nfft, d=sample_interval)
    )
    allowed = np.flatnonzero(np.abs(frequency_full) <= maximum_frequency_hz)
    frequency = frequency_full[allowed]
    spectrum = np.fft.fftshift(
        np.fft.fft(voltage, n=nfft, axis=0),
        axes=0,
    )[allowed] / len(times)
    inverse_noise = 1.0 / noise_power
    incoherent = np.sum(
        np.abs(spectrum) ** 2 * inverse_noise[None, :, None],
        axis=1,
    )
    feasible_by_range, altitude_by_range = grid.geometry_for_ranges(range_km)

    best_frequency = np.full(len(range_km), np.nan)
    best_velocity = np.full(len(range_km), np.nan)
    best_power_ratio = np.full(len(range_km), np.nan)
    candidate_rows = []
    denominator = float(np.sum(inverse_noise))

    for range_index in range(len(range_km)):
        feasible = feasible_by_range[range_index]
        if not np.any(feasible):
            continue
        direction_indices = np.flatnonzero(feasible)
        steering = grid.steering[direction_indices]
        preliminary = []
        best_preliminary_power = 0.0
        relative_threshold = 10.0 ** (
            AMBIGUITY_RELATIVE_POWER_DB / 10.0
        )
        doppler_indices = _top_doppler_indices(
            incoherent[:, range_index],
        )
        for doppler_index in doppler_indices:
            # Incoherent whitened power is an upper bound on every coherent
            # beam at this Doppler by Cauchy-Schwarz. Once this bound is below
            # the ambiguity threshold, no untested Doppler can be retained.
            upper_bound = float(incoherent[doppler_index, range_index])
            if (
                best_preliminary_power > 0.0
                and upper_bound < best_preliminary_power * relative_threshold
            ):
                break
            beam = np.einsum(
                "dc,c,c->d",
                steering,
                spectrum[doppler_index, :, range_index],
                inverse_noise,
                optimize=True,
            )
            beam_power = np.abs(beam) ** 2 / denominator
            local_in_feasible = _direction_local_maxima(
                grid,
                beam_power,
                feasible,
            )
            for grid_index in local_in_feasible:
                power = float(
                    beam_power[np.searchsorted(direction_indices, grid_index)]
                )
                best_preliminary_power = max(best_preliminary_power, power)
                preliminary.append(
                    (
                        power,
                        int(doppler_index),
                        int(grid_index),
                    )
                )
        if not preliminary:
            continue
        preliminary.sort(reverse=True)
        threshold = preliminary[0][0] * relative_threshold
        preliminary = [
            row for row in preliminary if row[0] >= threshold
        ][:MAX_AMBIGUITIES_PER_RANGE]

        refined = []
        for _, doppler_index, direction_index in preliminary:
            steering_vector = grid.steering[direction_index]
            neighbor_indices = np.asarray(
                [
                    max(0, doppler_index - 1),
                    doppler_index,
                    min(len(frequency) - 1, doppler_index + 1),
                ]
            )
            neighbor_beam = np.einsum(
                "c,kc,c->k",
                steering_vector,
                spectrum[neighbor_indices, :, range_index],
                inverse_noise,
                optimize=True,
            )
            fitted_frequency = _parabolic_frequency(
                frequency,
                np.abs(neighbor_beam) ** 2,
                doppler_index,
            )
            model = np.exp(2j * np.pi * fitted_frequency * times)
            amplitude = np.mean(
                np.conj(model)[:, None] * voltage[:, :, range_index],
                axis=0,
            )
            coherent_sum = np.sum(
                steering_vector * amplitude * inverse_noise
            )
            power_ratio = float(np.abs(coherent_sum) ** 2 / denominator)
            normalized_amplitude = amplitude / np.maximum(
                np.abs(amplitude), np.finfo(float).tiny
            )
            aoa_match = float(
                np.abs(np.mean(steering_vector * normalized_amplitude))
            )
            refined.append(
                (
                    power_ratio,
                    fitted_frequency,
                    direction_index,
                    aoa_match,
                )
            )
        refined.sort(reverse=True)
        maximum_power = refined[0][0]
        refined_threshold = maximum_power * 10.0 ** (
            AMBIGUITY_RELATIVE_POWER_DB / 10.0
        )
        refined = [
            row for row in refined if row[0] >= refined_threshold
        ][:MAX_AMBIGUITIES_PER_RANGE]
        best_power_ratio[range_index] = maximum_power
        best_frequency[range_index] = refined[0][1]
        best_velocity[range_index] = (
            -mc.wavelength * refined[0][1] / 2.0
        )

        for power_ratio, fitted_frequency, direction_index, aoa_match in refined:
            relative_db = 10.0 * np.log10(
                max(power_ratio, 1e-30) / max(maximum_power, 1e-30)
            )
            latitude, longitude, altitude_m = grid.geographic_position(
                range_km[range_index],
                direction_index,
            )
            candidate_rows.append(
                [
                    range_index,
                    range_km[range_index],
                    fitted_frequency,
                    -mc.wavelength * fitted_frequency / 2.0,
                    grid.u[direction_index],
                    grid.v[direction_index],
                    grid.w[direction_index],
                    latitude,
                    longitude,
                    altitude_m / 1000.0,
                    power_ratio,
                    relative_db,
                    aoa_match,
                ]
            )

    candidates = np.asarray(candidate_rows, dtype=np.float64)
    if candidates.size == 0:
        candidates = np.empty((0, len(CANDIDATE_COLUMNS)), dtype=np.float64)
    return best_frequency, best_velocity, best_power_ratio, candidates
