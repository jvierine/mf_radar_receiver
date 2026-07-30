#!/usr/bin/env python3
"""Plot ambiguity and provisional wind-fit diagnostics from realtime AoA fits."""

from __future__ import annotations

import os
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np


POWER_MIN_DB = 10.0
AOA_MATCH_MIN = 0.80
ALTITUDE_MIN_KM = 70.0
ALTITUDE_MAX_KM = 150.0
VELOCITY_LIMIT_MS = 150.0
WIND_KNOTS_KM = np.arange(70.0, 151.0, 10.0)
POWER_PENALTY_MS2_PER_DB = 25.0
SMOOTHNESS_WEIGHT = 20.0
HUBER_SCALE_MS = 20.0
RESIDUAL_ACCEPT_MS = 40.0

AMBIGUITY_PLOT = "latest_interferometry_ambiguities.png"
WIND_FIT_PLOT = "latest_interferometry_wind_fit.png"
DEBUG_DATA = "latest_interferometry_debug.h5"


def column_map(handle):
    return {
        name: index
        for index, name in enumerate(
            handle.attrs["aoa_candidate_columns"].split(",")
        )
    }


def linear_knot_basis(altitude_km):
    altitude_km = np.asarray(altitude_km)
    basis = np.zeros((len(altitude_km), len(WIND_KNOTS_KM)))
    upper = np.searchsorted(WIND_KNOTS_KM, altitude_km, side="right")
    upper = np.clip(upper, 1, len(WIND_KNOTS_KM) - 1)
    lower = upper - 1
    span = WIND_KNOTS_KM[upper] - WIND_KNOTS_KM[lower]
    upper_weight = (altitude_km - WIND_KNOTS_KM[lower]) / span
    basis[np.arange(len(basis)), lower] = 1.0 - upper_weight
    basis[np.arange(len(basis)), upper] = upper_weight
    return basis


def candidate_design(candidates, columns):
    basis = linear_knot_basis(candidates[:, columns["altitude_km"]])
    east = candidates[:, columns["east_direction_cosine"]]
    north = candidates[:, columns["north_direction_cosine"]]
    return np.hstack((basis * east[:, None], basis * north[:, None]))


def sorted_groups(candidates, columns):
    order = np.lexsort(
        (
            candidates[:, columns["range_index"]],
            candidates[:, columns["time_unix"]],
        )
    )
    candidates = candidates[order]
    keys = candidates[:, [columns["time_unix"], columns["range_index"]]]
    starts = np.r_[0, 1 + np.flatnonzero(np.any(np.diff(keys, axis=0), axis=1))]
    ends = np.r_[starts[1:], len(candidates)]
    return candidates, starts, ends


def select_per_cell(cost, starts, ends):
    return np.asarray(
        [start + np.argmin(cost[start:end]) for start, end in zip(starts, ends)],
        dtype=np.int64,
    )


def smoothness_rows(parameter_count):
    knot_count = len(WIND_KNOTS_KM)
    rows = []
    for offset in (0, knot_count):
        for index in range(knot_count - 2):
            row = np.zeros(parameter_count)
            row[offset + index : offset + index + 3] = (1.0, -2.0, 1.0)
            rows.append(row)
    return np.asarray(rows)


def ambiguity_resolving_wind_fit(candidates, columns, starts, ends):
    design = candidate_design(candidates, columns)
    velocity = candidates[:, columns["velocity_ms"]]
    relative_db = candidates[:, columns["relative_power_db"]]
    match = candidates[:, columns["aoa_match"]]
    power_db = 10.0 * np.log10(
        np.maximum(candidates[:, columns["beam_power_ratio"]], 1e-30)
    )

    beta = np.zeros(design.shape[1])
    selected = np.asarray([], dtype=np.int64)
    regularization = np.sqrt(SMOOTHNESS_WEIGHT) * smoothness_rows(
        design.shape[1]
    )
    for _ in range(12):
        predicted = design @ beta
        cost = (
            (velocity - predicted) ** 2
            + POWER_PENALTY_MS2_PER_DB * np.maximum(0.0, -relative_db)
            + 100.0 * (1.0 - match) ** 2
        )
        updated = select_per_cell(cost, starts, ends)
        x = design[updated]
        y = velocity[updated]
        base_weight = (
            np.clip(power_db[updated] - POWER_MIN_DB + 1.0, 1.0, 20.0)
            * match[updated] ** 2
        )
        residual = y - x @ beta
        robust_weight = np.minimum(
            1.0,
            HUBER_SCALE_MS / np.maximum(np.abs(residual), 1e-6),
        )
        weight = np.sqrt(base_weight * robust_weight)
        augmented_x = np.vstack((x * weight[:, None], regularization))
        augmented_y = np.r_[y * weight, np.zeros(len(regularization))]
        updated_beta = np.linalg.lstsq(
            augmented_x,
            augmented_y,
            rcond=None,
        )[0]
        if np.array_equal(updated, selected) and np.linalg.norm(
            updated_beta - beta
        ) < 0.05:
            beta = updated_beta
            selected = updated
            break
        beta = updated_beta
        selected = updated

    predicted = design[selected] @ beta
    residual = velocity[selected] - predicted
    accepted = np.abs(residual) <= RESIDUAL_ACCEPT_MS
    return beta, selected, predicted, residual, accepted


def style_axes(axes):
    for axis in np.ravel(axes):
        axis.set_facecolor("#050810")
        axis.tick_params(colors="#8fa1ba")
        axis.xaxis.label.set_color("#b7c5d9")
        axis.yaxis.label.set_color("#b7c5d9")
        axis.title.set_color("#edf4ff")
        for spine in axis.spines.values():
            spine.set_color("#314563")


def add_colorbar(figure, image, axis, label):
    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label(label, color="#b7c5d9")
    colorbar.ax.tick_params(colors="#8fa1ba")
    colorbar.outline.set_edgecolor("#314563")
    return colorbar


def save_figure(figure, destination):
    temporary = destination.with_suffix(".tmp.png")
    figure.savefig(
        temporary,
        dpi=160,
        facecolor=figure.get_facecolor(),
    )
    plt.close(figure)
    os.replace(temporary, destination)


def plot_ambiguities(candidates, columns, starts, ends, destination):
    counts = ends - starts
    representative = candidates[starts]
    sample_step = max(1, int(np.ceil(len(candidates) / 60_000)))
    sample = candidates[::sample_step]
    sample_range = sample[:, columns["range_km"]]
    sample_velocity = sample[:, columns["velocity_ms"]]
    east_km = sample_range * sample[:, columns["east_direction_cosine"]]
    north_km = sample_range * sample[:, columns["north_direction_cosine"]]

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(14, 10),
        facecolor="#070b14",
        constrained_layout=True,
    )
    style_axes(axes)
    time_values = [
        np.datetime64(int(value), "s").astype("datetime64[ms]").astype(object)
        for value in representative[:, columns["time_unix"]]
    ]
    count_plot = axes[0, 0].scatter(
        time_values,
        representative[:, columns["range_km"]],
        c=counts,
        s=5,
        cmap="viridis",
        vmin=1,
        vmax=12,
    )
    axes[0, 0].set_title("Retained joint Doppler–AoA candidates per cell")
    axes[0, 0].set_xlabel("Time (UTC)")
    axes[0, 0].set_ylabel("Round-trip range (km)")
    axes[0, 0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    add_colorbar(figure, count_plot, axes[0, 0], "Candidate count")

    position_plot = axes[0, 1].scatter(
        east_km,
        north_km,
        c=sample_velocity,
        s=2,
        alpha=0.35,
        cmap="seismic",
        vmin=-VELOCITY_LIMIT_MS,
        vmax=VELOCITY_LIMIT_MS,
        rasterized=True,
    )
    axes[0, 1].set_title("All retained ambiguous horizontal positions")
    axes[0, 1].set_xlabel("East (km)")
    axes[0, 1].set_ylabel("North (km)")
    axes[0, 1].set_aspect("equal", adjustable="box")
    add_colorbar(
        figure,
        position_plot,
        axes[0, 1],
        "Radial velocity (m/s)",
    )

    altitude_plot = axes[1, 0].scatter(
        sample_range,
        sample[:, columns["altitude_km"]],
        c=sample[:, columns["relative_power_db"]],
        s=2,
        alpha=0.4,
        cmap="magma",
        vmin=-6,
        vmax=0,
        rasterized=True,
    )
    axes[1, 0].set_title("WGS84 altitude allowed by each range")
    axes[1, 0].set_xlabel("Round-trip range (km)")
    axes[1, 0].set_ylabel("Altitude (km)")
    add_colorbar(figure, altitude_plot, axes[1, 0], "Relative power (dB)")

    axes[1, 1].hist(
        counts,
        bins=np.arange(0.5, 13.5, 1.0),
        color="#61d0ff",
        edgecolor="#070b14",
    )
    axes[1, 1].set_title("Ambiguity multiplicity")
    axes[1, 1].set_xlabel("Candidates retained per time–range cell")
    axes[1, 1].set_ylabel("Cells")
    figure.suptitle(
        "Three-dipole interferometry · ambiguity audit",
        color="#edf4ff",
        fontsize=19,
        weight="semibold",
    )
    save_figure(figure, destination)


def plot_wind_fit(
    candidates,
    columns,
    beta,
    selected,
    predicted,
    residual,
    accepted,
    destination,
):
    chosen = candidates[selected]
    accepted_rows = chosen[accepted]
    accepted_velocity = accepted_rows[:, columns["velocity_ms"]]
    accepted_prediction = predicted[accepted]
    accepted_residual = residual[accepted]
    accepted_range = accepted_rows[:, columns["range_km"]]
    east_km = (
        accepted_range
        * accepted_rows[:, columns["east_direction_cosine"]]
    )
    north_km = (
        accepted_range
        * accepted_rows[:, columns["north_direction_cosine"]]
    )
    knot_count = len(WIND_KNOTS_KM)

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(14, 10),
        facecolor="#070b14",
        constrained_layout=True,
    )
    style_axes(axes)
    position_plot = axes[0, 0].scatter(
        east_km,
        north_km,
        c=accepted_velocity,
        s=7,
        alpha=0.65,
        cmap="seismic",
        vmin=-VELOCITY_LIMIT_MS,
        vmax=VELOCITY_LIMIT_MS,
    )
    axes[0, 0].set_title("Automatically selected positions")
    axes[0, 0].set_xlabel("East (km)")
    axes[0, 0].set_ylabel("North (km)")
    axes[0, 0].set_aspect("equal", adjustable="box")
    add_colorbar(
        figure,
        position_plot,
        axes[0, 0],
        "Radial velocity (m/s)",
    )

    comparison_plot = axes[0, 1].scatter(
        accepted_prediction,
        accepted_velocity,
        c=accepted_rows[:, columns["altitude_km"]],
        s=7,
        alpha=0.6,
        cmap="viridis",
    )
    axes[0, 1].plot(
        [-VELOCITY_LIMIT_MS, VELOCITY_LIMIT_MS],
        [-VELOCITY_LIMIT_MS, VELOCITY_LIMIT_MS],
        color="#edf4ff",
        linewidth=1,
    )
    axes[0, 1].set_xlim(-VELOCITY_LIMIT_MS, VELOCITY_LIMIT_MS)
    axes[0, 1].set_ylim(-VELOCITY_LIMIT_MS, VELOCITY_LIMIT_MS)
    axes[0, 1].set_title("Observed versus wind-predicted radial velocity")
    axes[0, 1].set_xlabel("Predicted (m/s)")
    axes[0, 1].set_ylabel("Observed (m/s)")
    add_colorbar(
        figure,
        comparison_plot,
        axes[0, 1],
        "WGS84 altitude (km)",
    )

    axes[1, 0].plot(
        beta[:knot_count],
        WIND_KNOTS_KM,
        "o-",
        color="#ffb454",
        label="Zonal U",
    )
    axes[1, 0].plot(
        beta[knot_count:],
        WIND_KNOTS_KM,
        "o-",
        color="#61d0ff",
        label="Meridional V",
    )
    axes[1, 0].axvline(0, color="#8fa1ba", linewidth=0.8)
    axes[1, 0].set_title("Provisional altitude-smooth horizontal wind")
    axes[1, 0].set_xlabel("Wind velocity (m/s)")
    axes[1, 0].set_ylabel("WGS84 altitude (km)")
    axes[1, 0].legend(facecolor="#101827", labelcolor="#edf4ff")

    axes[1, 1].hist(
        accepted_residual,
        bins=np.linspace(-RESIDUAL_ACCEPT_MS, RESIDUAL_ACCEPT_MS, 41),
        color="#b786ff",
        edgecolor="#070b14",
    )
    rms = np.sqrt(np.mean(accepted_residual**2)) if len(accepted_residual) else np.nan
    axes[1, 1].set_title(
        f"Accepted radial residuals · N={len(accepted_rows)} · RMS={rms:.1f} m/s"
    )
    axes[1, 1].set_xlabel("Observed − predicted (m/s)")
    axes[1, 1].set_ylabel("Selections")
    figure.suptitle(
        "Diagnostic ambiguity-resolving zonal/meridional fit · not validated",
        color="#edf4ff",
        fontsize=19,
        weight="semibold",
    )
    save_figure(figure, destination)


def write_debug_data(
    destination,
    candidates,
    columns,
    beta,
    selected,
    predicted,
    residual,
    accepted,
):
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.attrs["method"] = (
            "alternating_discrete_ambiguity_selection_and_robust_"
            "altitude_smooth_horizontal_wind_fit"
        )
        handle.attrs["status"] = "diagnostic_not_validated"
        handle.attrs["candidate_columns"] = ",".join(columns)
        handle.attrs["power_min_db"] = POWER_MIN_DB
        handle.attrs["aoa_match_min"] = AOA_MATCH_MIN
        handle.attrs["residual_accept_ms"] = RESIDUAL_ACCEPT_MS
        handle.create_dataset("wind_altitude_km", data=WIND_KNOTS_KM)
        handle.create_dataset(
            "zonal_wind_ms",
            data=beta[: len(WIND_KNOTS_KM)],
        )
        handle.create_dataset(
            "meridional_wind_ms",
            data=beta[len(WIND_KNOTS_KM) :],
        )
        handle.create_dataset(
            "selected_candidates",
            data=candidates[selected],
            compression="gzip",
            compression_opts=1,
        )
        handle.create_dataset("predicted_radial_velocity_ms", data=predicted)
        handle.create_dataset("radial_velocity_residual_ms", data=residual)
        handle.create_dataset("accepted", data=accepted)
    os.replace(temporary, destination)


def main(data_file=None, plot_dir=None):
    data_file = Path(data_file or "/data2/plots/monitor/latest_doppler_15m.h5")
    plot_dir = Path(plot_dir or "/data2/plots/monitor")
    with h5py.File(data_file, "r") as handle:
        columns = column_map(handle)
        all_candidates = handle["aoa_candidates"][:]

    power_db = 10.0 * np.log10(
        np.maximum(all_candidates[:, columns["beam_power_ratio"]], 1e-30)
    )
    quality = (
        (power_db >= POWER_MIN_DB)
        & (all_candidates[:, columns["aoa_match"]] >= AOA_MATCH_MIN)
        & (all_candidates[:, columns["altitude_km"]] >= ALTITUDE_MIN_KM)
        & (all_candidates[:, columns["altitude_km"]] <= ALTITUDE_MAX_KM)
        & (
            np.abs(all_candidates[:, columns["velocity_ms"]])
            <= VELOCITY_LIMIT_MS
        )
    )
    audit_candidates, audit_starts, audit_ends = sorted_groups(
        all_candidates,
        columns,
    )
    candidates = all_candidates[quality]
    if not len(candidates):
        raise RuntimeError("no candidates pass the interferometry debug gates")
    candidates, starts, ends = sorted_groups(candidates, columns)
    beta, selected, predicted, residual, accepted = (
        ambiguity_resolving_wind_fit(candidates, columns, starts, ends)
    )
    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_ambiguities(
        audit_candidates,
        columns,
        audit_starts,
        audit_ends,
        plot_dir / AMBIGUITY_PLOT,
    )
    plot_wind_fit(
        candidates,
        columns,
        beta,
        selected,
        predicted,
        residual,
        accepted,
        plot_dir / WIND_FIT_PLOT,
    )
    write_debug_data(
        plot_dir / DEBUG_DATA,
        candidates,
        columns,
        beta,
        selected,
        predicted,
        residual,
        accepted,
    )
    print(
        f"Interferometry debug: {len(candidates)} candidates in {len(starts)}"
        f" cells; {np.count_nonzero(accepted)} provisional wind selections"
    )


if __name__ == "__main__":
    main()
