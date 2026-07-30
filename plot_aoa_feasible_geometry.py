#!/usr/bin/env python3
"""Plot WGS84 altitude limits for the MF-radar direction-cosine search."""

import argparse

import jcoord
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


RADAR_LAT = 69.58204
RADAR_LON = 19.22283
RADAR_ALT_M = 0.0
ALTITUDE_MIN_KM = 50.0
ALTITUDE_MAX_KM = 150.0
MINIMUM_ELEVATION_DEG = 20.0


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output")
    return parser.parse_args()


def wgs84_altitude_km(range_km, zenith_angle_deg, azimuth_deg=0.0):
    """Return geodetic altitude using jcoord's WGS84 conversion."""
    return (
        jcoord.az_el_r2geodetic(
            RADAR_LAT,
            RADAR_LON,
            RADAR_ALT_M,
            azimuth_deg,
            90.0 - zenith_angle_deg,
            range_km * 1000.0,
        )[2]
        / 1000.0
    )


def main():
    args = parse_args()
    ranges = np.linspace(50.0, 250.0, 401)
    zenith = np.linspace(0.0, 90.0 - MINIMUM_ELEVATION_DEG, 351)
    altitude = np.empty((len(zenith), len(ranges)))
    for range_index, range_km in enumerate(ranges):
        for angle_index, angle in enumerate(zenith):
            altitude[angle_index, range_index] = wgs84_altitude_km(
                range_km, angle
            )

    feasible = (
        (altitude >= ALTITUDE_MIN_KM)
        & (altitude <= ALTITUDE_MAX_KM)
    )
    figure, axis = plt.subplots(figsize=(10.5, 5.0), constrained_layout=True)
    axis.contourf(
        ranges,
        zenith,
        feasible.astype(float),
        levels=(-0.5, 0.5, 1.5),
        colors=("#eceff3", "#5b8cff"),
        alpha=0.95,
    )
    contours = axis.contour(
        ranges,
        zenith,
        altitude,
        levels=(50.0, 100.0, 150.0),
        colors=("#444444", "#ffffff", "#202020"),
        linewidths=(1.4, 1.2, 1.4),
    )
    axis.clabel(
        contours,
        fmt={50.0: "50 km", 100.0: "100 km", 150.0: "150 km"},
        inline=True,
        fontsize=9,
    )
    axis.set_xlim(50.0, 250.0)
    axis.set_ylim(0.0, 90.0 - MINIMUM_ELEVATION_DEG)
    axis.set_xlabel("Radar range coordinate (km)")
    axis.set_ylabel("Zenith angle (deg)")
    axis.set_title(
        "WGS84 AoA search region for 50–150 km meteor altitude"
    )
    axis.text(
        244,
        4,
        "blue: feasible",
        ha="right",
        va="bottom",
        color="#17243a",
    )
    axis.grid(color="#ffffff", linewidth=0.7, alpha=0.35)
    figure.savefig(args.output, dpi=180, facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    main()
