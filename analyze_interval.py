#!/usr/bin/env python3
"""Analyze a UTC interval with coherent RTIs and automatic wind detections."""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import json
import os
from pathlib import Path

import digital_rf as drf
import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import convolve1d

import mf_conf as mc
import calc_rti as crti
import plot_monitor_rti as monitor
import wind_estimates_3ch_2days as winds


DEFAULT_OUTPUT_ROOT = "/data2/products/interval_analysis"
WIND_BLOCK_S = 10 * 60
WIND_SEED_S = 10 * 60
POWER_THRESHOLD_DB = 6.0
REDUCTION_BLOCK_S = 2
REDUCTION_CHANNELS = ("ch1", "ch2", "ch3", "ch4")
_worker_reader = None
_worker_lpf = None
_worker_range_mask = None


def parse_utc(value: str) -> dt.datetime:
    parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    return parsed.astimezone(dt.timezone.utc)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("start", type=parse_utc, help="UTC interval start.")
    parser.add_argument("end", type=parse_utc, help="UTC interval end.")
    parser.add_argument("--raw-dir", default=mc.raw_dir)
    parser.add_argument("--metadata-dir", default=mc.xc_dir)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--cadence", type=int, default=60)
    parser.add_argument("--averages", type=int, default=5)
    parser.add_argument(
        "--materialize-metadata",
        action="store_true",
        help="Reduce the interval from raw voltage into an isolated metadata stream.",
    )
    parser.add_argument("--workers", type=int, default=4)
    return parser.parse_args()


def initialize_reduction_worker(raw_dir: str) -> None:
    global _worker_reader, _worker_lpf, _worker_range_mask
    _worker_reader = drf.DigitalRFReader(raw_dir)
    _worker_lpf = mc.fir_lowpass_hann(
        fc=20e3,
        fs=crti.FS_RAW,
        num_taps=crti.NUM_TAPS,
    ).astype(np.float32)
    _worker_range_mask = crti.range_mask(crti.range_vector())


def reduce_block(i0: int) -> tuple[int, dict]:
    if (
        _worker_reader is None
        or _worker_lpf is None
        or _worker_range_mask is None
    ):
        raise RuntimeError("reduction worker has no Digital RF reader")
    rtis = []
    rdis = []
    voltage_blocks = {
        channel: np.asarray(
            _worker_reader.read_vector_c81d(
                i0 + crti.OFFSET,
                REDUCTION_BLOCK_S * monitor.FS,
                channel,
            ),
            dtype=np.complex64,
        ).reshape(-1, crti.IPP)
        - np.complex64(mc.dc_offset)
        for channel in REDUCTION_CHANNELS
    }
    tx_phase = np.angle(
        np.mean(
            voltage_blocks["ch1"][:, : crti.TX_LEN],
            axis=1,
        )
    ).astype(np.float32)
    phase_correction = np.exp(-1j * tx_phase).astype(np.complex64)

    for channel in REDUCTION_CHANNELS:
        voltage = voltage_blocks[channel]
        voltage[:, : crti.GC] = 0.0
        voltage *= phase_correction[:, None]
        filtered = convolve1d(
            voltage,
            _worker_lpf,
            axis=1,
            mode="constant",
            cval=0.0,
            origin=-1,
        )
        decimated = filtered[:, :: crti.DEC]
        n_integrations = len(decimated) // crti.N_CI
        rti = decimated[: n_integrations * crti.N_CI].reshape(
            n_integrations,
            crti.N_CI,
            -1,
        ).sum(axis=1)
        rti = rti[:, _worker_range_mask].astype(np.complex64)
        rdi_full = np.fft.fftshift(
            np.fft.fft(rti, axis=0),
            axes=0,
        )
        rtis.append(rti)
        rdis.append(rdi_full.astype(np.complex64))

    tvec = (
        np.arange(n_integrations)
        * crti.N_CI
        * crti.IPP
        / crti.FS_RAW
    )
    rvec = crti.range_vector()[_worker_range_mask]
    prf_integrated = crti.FS_RAW / crti.IPP / crti.N_CI
    fvec = np.fft.fftshift(
        np.fft.fftfreq(n_integrations, d=1.0 / prf_integrated)
    )
    return i0, {
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


def materialize_interval_metadata(
    raw_dir: str,
    metadata_dir: Path,
    start_unix: int,
    end_unix: int,
    workers: int,
) -> None:
    metadata_dir.mkdir(parents=True, exist_ok=True)
    writer = drf.DigitalMetadataWriter(
        str(metadata_dir),
        3600,
        REDUCTION_BLOCK_S,
        monitor.FS,
        1,
        "xc",
    )
    starts = np.arange(
        start_unix * monitor.FS,
        end_unix * monitor.FS,
        REDUCTION_BLOCK_S * monitor.FS,
        dtype=np.int64,
    )
    batch_size = max(1, workers * 10)
    completed = 0
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=workers,
        initializer=initialize_reduction_worker,
        initargs=(raw_dir,),
    ) as executor:
        for batch_start in range(0, len(starts), batch_size):
            batch = starts[batch_start : batch_start + batch_size]
            for i0, value in executor.map(reduce_block, map(int, batch)):
                writer.write(i0, value)
                completed += 1
            print(
                f"Reduced {completed}/{len(starts)} two-second blocks",
                flush=True,
            )


def band_statistics(
    power_ratio_db: np.ndarray,
    range_km: np.ndarray,
    lower_km: float,
    upper_km: float,
) -> dict:
    mask = (range_km >= lower_km) & (range_km < upper_km)
    values = power_ratio_db[:, mask]
    profile_peaks = np.nanmax(values, axis=1)
    persistence = np.nanmean(values >= POWER_THRESHOLD_DB, axis=0)
    band_ranges = range_km[mask]
    persistent_ranges = band_ranges[persistence >= 0.5]
    return {
        "range_km": [lower_km, upper_km],
        "threshold_db": POWER_THRESHOLD_DB,
        "median_profile_peak_db": float(np.nanmedian(profile_peaks)),
        "maximum_db": float(np.nanmax(values)),
        "profiles_with_peak_above_threshold_fraction": float(
            np.nanmean(profile_peaks >= POWER_THRESHOLD_DB)
        ),
        "cells_above_threshold_fraction": float(
            np.nanmean(values >= POWER_THRESHOLD_DB)
        ),
        "maximum_gate_persistence_fraction": float(np.nanmax(persistence)),
        "persistent_gate_ranges_km": persistent_ranges.round(2).tolist(),
    }


def automatic_detections(
    metadata_reader: drf.DigitalMetadataReader,
    start_unix: float,
    end_unix: float,
) -> tuple[np.ndarray, np.ndarray]:
    phasecal = winds.load_phasecal()
    coordinates = winds.antenna_coordinates_ecef()
    position_differences = [
        coordinates[0, :] - coordinates[2, :],
        coordinates[0, :] - coordinates[3, :],
        coordinates[2, :] - coordinates[3, :],
    ]

    detections = []
    wind_rows = []
    prior_profile = None
    processing_start = start_unix - WIND_SEED_S

    for block_start in np.arange(processing_start, end_unix, WIND_BLOCK_S):
        block_end = min(block_start + WIND_BLOCK_S, end_unix)
        block_detections = winds.extract_detections_for_interval(
            metadata_reader,
            phasecal,
            position_differences,
            block_start,
            block_end,
            prior_wind_profile=prior_profile,
        )
        block_winds = (
            winds.estimate_wind_profile(block_detections, block_start)
            if len(block_detections)
            else winds.empty_wind_block(block_start)
        )
        wind_rows.append(block_winds)
        prior_profile = winds.latest_accepted_wind_profile(
            np.vstack(wind_rows)
        )
        if block_start >= start_unix and len(block_detections):
            detections.append(block_detections)

    selected = (
        np.vstack(detections)
        if detections
        else np.empty((0, 21), dtype=np.float64)
    )
    selected = winds.deduplicate_detections(selected)
    requested_winds = np.vstack(wind_rows)[
        np.vstack(wind_rows)[:, 0] >= start_unix
    ]
    return selected, requested_winds


def selection_funnel(
    metadata_reader: drf.DigitalMetadataReader,
    start_unix: float,
    end_unix: float,
) -> dict:
    """Count candidates surviving each sequential automatic selection gate."""
    phasecal = winds.load_phasecal()
    counts = {
        "two_second_records": 0,
        "range_doppler_cells": 0,
        "after_snr": 0,
        "after_background_rejection": 0,
        "after_coherence": 0,
        "after_phase_closure": 0,
        "locator_pixels_after_one_per_range_gate": 0,
        "range_gates_tested_for_background": 0,
        "range_gates_rejected_as_broad_background": 0,
    }
    t0_us = int(start_unix * 1e6)
    t1_us = int(end_unix * 1e6)

    for read_start in np.arange(t0_us, t1_us, int(winds.READ_DT * 1e6)):
        read_end = min(read_start + int(winds.READ_DT * 1e6), t1_us)
        for record in metadata_reader.read(read_start, read_end).values():
            if not all(
                name in record
                for name in ("rdi1", "rdi3", "rdi4", "rvec", "fvec")
            ):
                continue
            counts["two_second_records"] += 1
            rdi1 = record["rdi1"] * np.exp(1j * phasecal[0])
            rdi3 = record["rdi3"] * np.exp(1j * phasecal[2])
            rdi4 = record["rdi4"] * np.exp(1j * phasecal[3])
            rvec = record["rvec"]
            fvec = record["fvec"]

            try:
                ri0 = np.where(rvec > winds.NOISER0)[0][0]
                ri1 = np.where(rvec > winds.NOISER1)[0][0]
                fi0 = np.where(fvec < -winds.NOISE_FMIN)[0][-1]
                fi1 = np.where(fvec > winds.NOISE_FMIN)[0][0]
            except IndexError:
                continue
            noise_power = 0.5 * (
                np.mean(np.abs(rdi1[:fi0, ri0:ri1]) ** 2)
                + np.mean(np.abs(rdi1[fi1:-1, ri0:ri1]) ** 2)
            )
            if noise_power <= 0 or not np.isfinite(noise_power):
                continue

            snr = (np.abs(rdi1) ** 2 - noise_power) / noise_power
            rr, ff = np.meshgrid(rvec, fvec)
            geometry = (
                (rr > winds.RANGE_MIN)
                & (rr < winds.RANGE_MAX)
                & (np.abs(ff) <= winds.MAX_FIT_DOPPLER_HZ)
            )
            snr_pass = geometry & (snr > winds.SNR_THRESH)

            strong_fraction = np.mean(snr > winds.SNR_THRESH, axis=0)
            good_gate = strong_fraction < winds.BG_STRONG_FRAC_MAX
            background_pass = snr_pass & good_gate[None, :]

            cross_spectra, coherence = winds.three_dipole_coherence(
                rdi1,
                rdi3,
                rdi4,
            )
            coherence_pass = background_pass & (
                coherence >= winds.COHERENCE_MIN
            )
            closure = winds.phase_closure(*cross_spectra)
            closure_pass = coherence_pass & (
                closure <= winds.PHASE_CLOSURE_MAX
            )

            counts["range_doppler_cells"] += int(np.count_nonzero(geometry))
            counts["after_snr"] += int(np.count_nonzero(snr_pass))
            counts["after_background_rejection"] += int(
                np.count_nonzero(background_pass)
            )
            counts["after_coherence"] += int(
                np.count_nonzero(coherence_pass)
            )
            counts["after_phase_closure"] += int(
                np.count_nonzero(closure_pass)
            )
            tested_gate = (rvec > winds.RANGE_MIN) & (rvec < winds.RANGE_MAX)
            counts["range_gates_tested_for_background"] += int(
                np.count_nonzero(tested_gate)
            )
            counts["range_gates_rejected_as_broad_background"] += int(
                np.count_nonzero(tested_gate & ~good_gate)
            )
            candidate_indices = np.where(closure_pass)
            counts["locator_pixels_after_one_per_range_gate"] += int(
                len(np.unique(candidate_indices[1]))
                if len(candidate_indices[0])
                else 0
            )

    return counts


def detection_statistics(detections: np.ndarray) -> dict:
    if len(detections) == 0:
        return {"count": 0}

    snr_db = 10.0 * np.log10(np.maximum(detections[:, 5], 1e-20))
    result = {
        "count": int(len(detections)),
        "unique_two_second_times": int(len(np.unique(detections[:, 0]))),
        "time_span_seconds": float(
            np.max(detections[:, 0]) - np.min(detections[:, 0])
        ),
        "height_km": {
            "minimum": float(np.min(detections[:, 3])),
            "median": float(np.median(detections[:, 3])),
            "maximum": float(np.max(detections[:, 3])),
        },
        "beamformed_snr_db": {
            "minimum": float(np.min(snr_db)),
            "median": float(np.median(snr_db)),
            "maximum": float(np.max(snr_db)),
        },
        "coherence": {
            "minimum": float(np.min(detections[:, 11])),
            "median": float(np.median(detections[:, 11])),
            "maximum": float(np.max(detections[:, 11])),
        },
        "radial_velocity_ms": {
            "minimum": float(np.min(detections[:, 4])),
            "median": float(np.median(detections[:, 4])),
            "maximum": float(np.max(detections[:, 4])),
        },
        "east_km": [
            float(np.min(detections[:, 1])),
            float(np.max(detections[:, 1])),
        ],
        "north_km": [
            float(np.min(detections[:, 2])),
            float(np.max(detections[:, 2])),
        ],
    }
    result["altitude_cut_counts"] = {
        f"{altitude:.0f}_plus_minus_2_km": int(
            np.count_nonzero(np.abs(detections[:, 3] - altitude) <= 2.0)
        )
        for altitude in winds.ALTITUDE_CUT_CENTERS_KM
    }
    return result


def wind_statistics(wind_rows: np.ndarray) -> dict:
    finite = np.isfinite(wind_rows[:, 2]) & np.isfinite(wind_rows[:, 3])
    return {
        "rows": int(len(wind_rows)),
        "finite_rows": int(np.count_nonzero(finite)),
        "finite_fraction": float(np.mean(finite)) if len(finite) else 0.0,
        "zonal_ms": (
            [
                float(np.min(wind_rows[finite, 2])),
                float(np.median(wind_rows[finite, 2])),
                float(np.max(wind_rows[finite, 2])),
            ]
            if np.any(finite)
            else None
        ),
        "meridional_ms": (
            [
                float(np.min(wind_rows[finite, 3])),
                float(np.median(wind_rows[finite, 3])),
                float(np.max(wind_rows[finite, 3])),
            ]
            if np.any(finite)
            else None
        ),
    }


def image_extent(
    times: np.ndarray,
    range_km: np.ndarray,
) -> tuple[float, float, float, float]:
    date_numbers = mdates.date2num(times.astype("datetime64[s]"))
    half_time_step = 0.5 * float(np.median(np.diff(date_numbers)))
    half_range_step = 0.5 * float(np.median(np.diff(range_km)))
    return (
        date_numbers[0] - half_time_step,
        date_numbers[-1] + half_time_step,
        range_km[0] - half_range_step,
        range_km[-1] + half_range_step,
    )


def plot_analysis(
    times: np.ndarray,
    range_km: np.ndarray,
    power_ratio_db: np.ndarray,
    detections: np.ndarray,
    start: dt.datetime,
    end: dt.datetime,
    output_path: Path,
) -> None:
    figure, axes = plt.subplots(
        3,
        1,
        figsize=(15, 13),
        facecolor="#070b14",
        constrained_layout=True,
    )
    for axis in axes:
        axis.set_facecolor("#050810")
        axis.tick_params(colors="#8fa1ba")
        axis.xaxis.label.set_color("#b7c5d9")
        axis.yaxis.label.set_color("#b7c5d9")
        axis.title.set_color("#edf4ff")
        for spine in axis.spines.values():
            spine.set_color("#314563")

    for axis, upper_range, limits, title in (
        (axes[0], 1500.0, (-3.0, 35.0), "Full-range coherent RTI"),
        (axes[1], 200.0, (-10.0, 20.0), "Mesospheric coherent RTI"),
    ):
        mask = (range_km >= 0) & (range_km <= upper_range)
        image = axis.imshow(
            power_ratio_db[:, mask].T,
            extent=image_extent(times, range_km[mask]),
            origin="lower",
            aspect="auto",
            interpolation="nearest",
            cmap="plasma",
            vmin=limits[0],
            vmax=limits[1],
        )
        axis.xaxis_date()
        axis.xaxis.set_major_formatter(
            mdates.DateFormatter("%H:%M", tz=dt.timezone.utc)
        )
        axis.set_ylim(0, upper_range)
        axis.set_ylabel("One-way range (km)")
        axis.set_title(title)
        colorbar = figure.colorbar(image, ax=axis, pad=0.01)
        colorbar.set_label("Power / 30–50 km background (dB)")

    if len(detections):
        axes[1].scatter(
            [
                dt.datetime.fromtimestamp(value, tz=dt.timezone.utc)
                for value in detections[:, 0]
            ],
            detections[:, 9],
            s=8,
            facecolors="none",
            edgecolors="#ffffff",
            linewidths=0.45,
            alpha=0.8,
        )

        velocity_limit = float(
            np.clip(np.percentile(np.abs(detections[:, 4]), 98), 20, 300)
        )
        positions = axes[2].scatter(
            detections[:, 1],
            detections[:, 2],
            c=detections[:, 4],
            s=10,
            cmap="seismic",
            vmin=-velocity_limit,
            vmax=velocity_limit,
            linewidths=0,
        )
        figure.colorbar(positions, ax=axes[2], pad=0.01).set_label(
            "Radial velocity (m/s)"
        )
    else:
        axes[2].text(
            0.5,
            0.5,
            "No automatic wind detections",
            transform=axes[2].transAxes,
            ha="center",
            va="center",
            color="#edf4ff",
        )

    axes[2].scatter([0], [0], marker="+", s=90, color="#edf4ff")
    axes[2].set_xlim(-200, 200)
    axes[2].set_ylim(-200, 200)
    axes[2].set_aspect("equal", adjustable="box")
    axes[2].set_xlabel("East–west position (km)")
    axes[2].set_ylabel("North–south position (km)")
    axes[2].set_title("Automatically accepted echo positions")
    axes[1].set_xlabel("Time (UTC)")
    axes[0].set_xlabel("Time (UTC)")
    figure.suptitle(
        f"Ramfjordmoen MF radar interval analysis\n"
        f"{start:%Y-%m-%d %H:%M}–{end:%H:%M} UTC",
        color="#edf4ff",
        fontsize=17,
    )
    temporary = output_path.with_suffix(".tmp.png")
    figure.savefig(temporary, dpi=160, facecolor=figure.get_facecolor())
    plt.close(figure)
    os.replace(temporary, output_path)


def main() -> None:
    args = parse_args()
    if args.end <= args.start:
        raise ValueError("end must be after start")

    start_unix = int(args.start.timestamp())
    end_unix = int(args.end.timestamp())
    times = np.arange(start_unix, end_unix, args.cadence, dtype=np.int64)
    if len(times) < 2:
        raise ValueError("interval must contain at least two cadence bins")

    tag = f"{args.start:%Y%m%dT%H%M}_{args.end:%H%M}"
    output_dir = Path(args.output_root) / tag
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_reader = drf.DigitalRFReader(args.raw_dir)
    power = np.asarray(
        [
            monitor.process_time_bin(
                raw_reader,
                int(timestamp),
                list(monitor.DEFAULT_CHANNELS),
                args.averages,
            )
            for timestamp in times
        ],
        dtype=np.float32,
    )
    range_km = monitor.range_vector_km(power.shape[1])
    background_mask = (
        (range_km >= monitor.THIRTY_MINUTE_NOISE_RANGE_KM[0])
        & (range_km <= monitor.THIRTY_MINUTE_NOISE_RANGE_KM[1])
    )
    background_power = float(np.mean(power[:, background_mask]))
    power_ratio_db = 10.0 * np.log10(
        np.maximum(power / background_power, 1e-20)
    )

    if args.materialize_metadata:
        interval_metadata_dir = output_dir / "xc"
        materialize_interval_metadata(
            args.raw_dir,
            interval_metadata_dir,
            start_unix - WIND_SEED_S,
            end_unix,
            args.workers,
        )
        metadata_dir = str(interval_metadata_dir)
    else:
        metadata_dir = args.metadata_dir

    metadata_reader = drf.DigitalMetadataReader(metadata_dir)
    detections, wind_rows = automatic_detections(
        metadata_reader,
        start_unix,
        end_unix,
    )

    plot_path = output_dir / "interval_analysis.png"
    plot_analysis(
        times,
        range_km,
        power_ratio_db,
        detections,
        args.start,
        args.end,
        plot_path,
    )
    np.save(output_dir / "detections.npy", detections)
    np.save(output_dir / "winds.npy", wind_rows)

    metrics = {
        "start": args.start.isoformat().replace("+00:00", "Z"),
        "end": args.end.isoformat().replace("+00:00", "Z"),
        "background_method": "mean_power_over_complete_interval_and_30_to_50_km",
        "background_power": background_power,
        "mesospheric_slant_range": band_statistics(
            power_ratio_db,
            range_km,
            60.0,
            200.0,
        ),
        "ionospheric_slant_range": band_statistics(
            power_ratio_db,
            range_km,
            200.0,
            1500.0,
        ),
        "selection_funnel": selection_funnel(
            metadata_reader,
            start_unix,
            end_unix,
        ),
        "automatic_detections": detection_statistics(detections),
        "wind_estimates": wind_statistics(wind_rows),
    }
    (output_dir / "metrics.json").write_text(
        json.dumps(metrics, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(metrics, indent=2, allow_nan=False))
    print(f"Plot: {plot_path}")


if __name__ == "__main__":
    main()
