#!/usr/bin/env python3
"""Estimate receiver-channel phase corrections from vertical radar echoes.

The direct transmit path is never used.  Candidate calibration samples are
high-SNR, high-coherence echoes at mesospheric ranges and in Doppler bins
closest to zero.  The ensemble mean arrival direction is assumed to be
vertical, so its expected geometric phase is zero on every baseline.
"""

import argparse
import datetime as dt
import os
import shutil
import tempfile

import h5py
import numpy as np


CHANNELS = (1, 2, 3, 4)
DIPOLE_INDICES = (0, 2, 3)
CALIBRATED_INDICES = DIPOLE_INDICES
REFERENCE_CHANNEL = 1
FIT_DURATION_S = 10
SOURCE_RECORD_S = 2


def circular_statistics(weighted_sum, weight_sum):
    """Return circular mean, resultant length, and circular standard deviation."""
    mean = np.angle(weighted_sum)
    resultant = np.clip(np.abs(weighted_sum) / np.maximum(weight_sum, 1e-30), 0, 1)
    circular_std = np.sqrt(np.maximum(0.0, -2.0 * np.log(np.maximum(resultant, 1e-12))))
    return mean, resultant, circular_std


def parse_time(value):
    """Parse Unix seconds or an ISO-8601 UTC time."""
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
        if parsed.tzinfo is None:
            parsed = parsed.replace(tzinfo=dt.timezone.utc)
        return parsed.timestamp()


def iter_fit_groups(reader, start_seconds, end_seconds):
    """Yield complete, non-overlapping groups spanning exactly ten seconds."""
    first_start = int(np.ceil(start_seconds / FIT_DURATION_S) * FIT_DURATION_S)
    final_start = int(np.floor(end_seconds / FIT_DURATION_S) * FIT_DURATION_S)
    for group_start in range(first_start, final_start, FIT_DURATION_S):
        keys = [
            int((group_start + offset) * 1e6)
            for offset in range(0, FIT_DURATION_S, SOURCE_RECORD_S)
        ]
        metadata = reader.read(
            keys[0],
            int((group_start + FIT_DURATION_S) * 1e6) - 1,
        )
        if not all(key in metadata for key in keys):
            continue
        records = [metadata[key] for key in keys]
        yield group_start + FIT_DURATION_S / 2.0, keys, records


def fit_group_phasors(
    keys,
    records,
    range_limits,
    doppler_max,
    snr_min_db,
    coherence_min,
):
    """Return fitted-amplitude phase phasors for one complete ten-second group."""
    from plot_dense_doppler import fit_common_sinusoid_fft
    import wind_estimates_3ch_2days as winds

    fields = [f"rti{channel}" for channel in CHANNELS]
    if not all(
        field in record
        for record in records
        for field in (*fields, "rvec", "tvec")
    ):
        return None

    rvec = np.asarray(records[0]["rvec"], dtype=np.float64)
    if not all(
        np.array_equal(rvec, np.asarray(record["rvec"], dtype=np.float64))
        for record in records[1:]
    ):
        return None
    range_mask = (rvec >= range_limits[0]) & (rvec <= range_limits[1])
    if not np.any(range_mask):
        return None

    fit_times = np.concatenate(
        [
            key / 1e6 + np.asarray(record["tvec"], dtype=np.float64)
            for key, record in zip(keys, records)
        ]
    )
    fit_times -= fit_times[0]
    channel_voltage = np.asarray(
        [
            np.concatenate(
                [
                    np.asarray(record[field])[:, range_mask]
                    for record in records
                ],
                axis=0,
            )
            for field in fields
        ],
        dtype=np.complex128,
    ).transpose(1, 0, 2)
    frequency, _, dipole_amplitude = fit_common_sinusoid_fft(
        fit_times,
        channel_voltage[:, np.asarray(DIPOLE_INDICES), :],
        winds.MAX_FIT_DOPPLER_HZ,
    )

    model = np.exp(2j * np.pi * fit_times[:, None] * frequency[None, :])
    amplitude = np.mean(
        np.conj(model)[:, None, :] * channel_voltage,
        axis=0,
    )
    amplitude[np.asarray(DIPOLE_INDICES)] = dipole_amplitude
    residual = channel_voltage - model[:, None, :] * amplitude[None, :, :]
    per_channel_snr_db = 10.0 * np.log10(
        np.maximum(
            np.abs(amplitude) ** 2
            / np.maximum(np.mean(np.abs(residual) ** 2, axis=0), 1e-30),
            1e-30,
        )
    )
    dipole_snr_db = np.min(
        per_channel_snr_db[np.asarray(DIPOLE_INDICES)],
        axis=0,
    )

    # Five independent two-second fitted-amplitude looks measure whether the
    # inter-channel phase remains coherent throughout the ten-second fit.
    samples_per_record = len(np.asarray(records[0]["tvec"]))
    segment_amplitude = np.empty(
        (len(records), len(CHANNELS), int(np.count_nonzero(range_mask))),
        dtype=np.complex128,
    )
    for segment_index in range(len(records)):
        sample_slice = slice(
            segment_index * samples_per_record,
            (segment_index + 1) * samples_per_record,
        )
        segment_amplitude[segment_index] = np.mean(
            np.conj(model[sample_slice])[:, None, :]
            * channel_voltage[sample_slice],
            axis=0,
        )

    pair_cross = {}
    pair_coherence = {}
    tiny = np.finfo(np.float32).tiny
    for first in range(len(CHANNELS)):
        for second in range(first + 1, len(CHANNELS)):
            cross = np.mean(
                segment_amplitude[:, first]
                * np.conj(segment_amplitude[:, second]),
                axis=0,
            )
            denominator = np.sqrt(
                np.maximum(
                    np.mean(np.abs(segment_amplitude[:, first]) ** 2, axis=0)
                    * np.mean(np.abs(segment_amplitude[:, second]) ** 2, axis=0),
                    tiny,
                )
            )
            pair_cross[(first, second)] = cross
            pair_coherence[(first, second)] = np.clip(
                np.abs(cross) / denominator, 0.0, 1.0
            )

    dipole_coherence = np.mean(
        [
            pair_coherence[(0, 2)],
            pair_coherence[(0, 3)],
            pair_coherence[(2, 3)],
        ],
        axis=0,
    )
    candidate = (
        (np.abs(frequency) <= doppler_max)
        & (dipole_snr_db >= snr_min_db)
        & (dipole_coherence >= coherence_min)
    )
    if not np.any(candidate):
        return None

    phasors = np.ones((len(CHANNELS), int(np.count_nonzero(candidate))), complex)
    coherences = np.ones_like(phasors.real)
    for channel_index in range(1, len(CHANNELS)):
        cross = amplitude[0, candidate] * np.conj(
            amplitude[channel_index, candidate]
        )
        phasors[channel_index] = cross / np.maximum(np.abs(cross), tiny)
        coherences[channel_index] = pair_coherence[(0, channel_index)][candidate]

    # Equalize the influence of very strong echoes.  Coherence remains a useful
    # reliability weight, while amplitude is used only for the SNR threshold.
    weights = coherences**2
    weights[0] = dipole_coherence[candidate] ** 2
    return phasors, weights


def write_calibration(path, result):
    """Atomically write a backward-compatible HDF5 calibration product."""
    directory = os.path.abspath(os.path.dirname(path) or ".")
    os.makedirs(directory, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=".phasecal-", suffix=".h5", dir=directory
    )
    os.close(fd)
    try:
        with h5py.File(temporary, "w") as output:
            for name, values in result["datasets"].items():
                output.create_dataset(name, data=values)
            for name, value in result["attrs"].items():
                output.attrs[name] = value
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def install_calibration(candidate_path, active_path, history_dir):
    """Archive the active HDF5 product, then atomically install the candidate."""
    os.makedirs(history_dir, exist_ok=True)
    if os.path.exists(active_path):
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        archive_path = os.path.join(history_dir, f"phasecal_before_{stamp}.h5")
        shutil.copy2(active_path, archive_path)
    directory = os.path.abspath(os.path.dirname(active_path) or ".")
    fd, temporary = tempfile.mkstemp(prefix=".phasecal-install-", dir=directory)
    os.close(fd)
    try:
        shutil.copyfile(candidate_path, temporary)
        os.replace(temporary, active_path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-dir", default="/data2/metadata/xc")
    parser.add_argument("--hours", type=float, default=24.0)
    parser.add_argument("--end", help="UTC ISO time or Unix seconds; default newest")
    parser.add_argument("--range-km", nargs=2, type=float, default=(70.0, 150.0))
    parser.add_argument("--doppler-max-hz", type=float, default=0.25)
    parser.add_argument("--snr-min-db", type=float, default=10.0)
    parser.add_argument("--coherence-min", type=float, default=0.80)
    parser.add_argument("--stability-block-minutes", type=int, default=30)
    parser.add_argument("--minimum-cells", type=int, default=1000)
    parser.add_argument("--minimum-resultant", type=float, default=0.50)
    parser.add_argument("--maximum-block-std-rad", type=float, default=0.35)
    parser.add_argument("--output", default="phasecal_candidate.h5")
    parser.add_argument(
        "--install",
        action="store_true",
        help="install the accepted candidate as the active phasecal.h5",
    )
    parser.add_argument("--active-path", default="phasecal.h5")
    parser.add_argument(
        "--history-dir", default="/data2/products/calibration/phasecal_history"
    )
    return parser


def main():
    args = build_parser().parse_args()
    import digital_rf as drf

    reader = drf.DigitalMetadataReader(args.metadata_dir)
    bounds = reader.get_bounds()
    end_seconds = parse_time(args.end)
    if end_seconds is None:
        end_seconds = bounds[1] / 1e6
    start_seconds = max(bounds[0] / 1e6, end_seconds - args.hours * 3600.0)
    if start_seconds >= end_seconds:
        raise ValueError("calibration interval contains no metadata")

    block_seconds = args.stability_block_minutes * 60
    block_edges = np.arange(start_seconds, end_seconds + block_seconds, block_seconds)
    if block_edges[-1] < end_seconds:
        block_edges = np.append(block_edges, end_seconds)
    block_sum = np.zeros((len(block_edges) - 1, len(CHANNELS)), complex)
    block_weight = np.zeros((len(block_edges) - 1, len(CHANNELS)), float)
    block_count = np.zeros(len(block_edges) - 1, np.int64)

    records_used = 0
    for center_time, keys, records in iter_fit_groups(
        reader, start_seconds, end_seconds
    ):
        selected = fit_group_phasors(
            keys,
            records,
            tuple(args.range_km),
            args.doppler_max_hz,
            args.snr_min_db,
            args.coherence_min,
        )
        if selected is None:
            continue
        phasors, weights = selected
        block_index = min(
            int((center_time - start_seconds) // block_seconds),
            len(block_count) - 1,
        )
        block_sum[block_index] += np.sum(phasors * weights, axis=1)
        block_weight[block_index] += np.sum(weights, axis=1)
        block_count[block_index] += phasors.shape[1]
        records_used += 1

    total_sum = np.sum(block_sum, axis=0)
    total_weight = np.sum(block_weight, axis=0)
    phasecal, resultant, circular_std = circular_statistics(total_sum, total_weight)
    phasecal[0] = 0.0
    calibrated_channel_mask = np.zeros(len(CHANNELS), dtype=np.bool_)
    calibrated_channel_mask[np.asarray(CALIBRATED_INDICES)] = True
    loop_phase_retained = False
    if os.path.exists(args.active_path):
        with h5py.File(args.active_path, "r") as active:
            active_phasecal = np.asarray(active["phasecal"][()])
        if active_phasecal.shape == phasecal.shape:
            phasecal[1] = active_phasecal[1]
            loop_phase_retained = True

    valid_blocks = block_count > 0
    block_phase = np.full_like(block_weight, np.nan)
    block_resultant = np.full_like(block_weight, np.nan)
    if np.any(valid_blocks):
        stats = circular_statistics(
            block_sum[valid_blocks], block_weight[valid_blocks]
        )
        block_phase[valid_blocks] = stats[0]
        block_resultant[valid_blocks] = stats[1]
    block_phase[:, 0] = 0.0

    phase_difference = np.angle(
        np.exp(1j * (block_phase[valid_blocks] - phasecal[None, :]))
    )
    block_std = np.nanstd(phase_difference, axis=0)
    block_std[0] = 0.0
    total_cells = int(np.sum(block_count))
    accepted = bool(
        total_cells >= args.minimum_cells
        and np.all(
            resultant[np.asarray(CALIBRATED_INDICES[1:])]
            >= args.minimum_resultant
        )
        and np.all(
            block_std[np.asarray(CALIBRATED_INDICES[1:])]
            <= args.maximum_block_std_rad
        )
    )

    now = dt.datetime.now(dt.timezone.utc)
    result = {
        "datasets": {
            "phasecal": phasecal,
            "calibrated_channel_mask": calibrated_channel_mask,
            "circular_resultant": resultant,
            "circular_std_rad": circular_std,
            "block_start_unix": block_edges[:-1],
            "block_phasecal_rad": block_phase,
            "block_resultant": block_resultant,
            "block_selected_cells": block_count,
            "block_phase_std_rad": block_std,
        },
        "attrs": {
            "created_utc": now.isoformat(),
            "method": "vertical_echo_joint_10s_sinusoid_fit",
            "direct_transmit_path_used": False,
            "fit_duration_seconds": FIT_DURATION_S,
            "source_record_seconds": SOURCE_RECORD_S,
            "reference_channel": REFERENCE_CHANNEL,
            "channels": np.asarray(CHANNELS),
            "loop_phase_retained_from_previous_calibration": loop_phase_retained,
            "interval_start_unix": start_seconds,
            "interval_end_unix": end_seconds,
            "range_min_km": args.range_km[0],
            "range_max_km": args.range_km[1],
            "doppler_abs_max_hz": args.doppler_max_hz,
            "snr_min_db": args.snr_min_db,
            "coherence_min": args.coherence_min,
            "selected_cells": total_cells,
            "records_used": records_used,
            "accepted": accepted,
            "acceptance_minimum_cells": args.minimum_cells,
            "acceptance_minimum_resultant": args.minimum_resultant,
            "acceptance_maximum_block_std_rad": args.maximum_block_std_rad,
            "assumption": (
                "The ensemble of high-SNR, high-coherence 10-second sinusoid "
                "fits closest to zero Doppler is centered on the vertical direction."
            ),
            "loop_status": (
                "diagnostic_only_not_calibrated; loop phase is not stable "
                "under the vertical-echo assumption"
            ),
        },
    }
    write_calibration(args.output, result)

    print(f"interval UTC: {dt.datetime.fromtimestamp(start_seconds, dt.timezone.utc)}")
    print(f"              {dt.datetime.fromtimestamp(end_seconds, dt.timezone.utc)}")
    print(f"records used: {records_used}")
    print(f"selected cells: {total_cells}")
    print(f"phasecal rad: {np.array2string(phasecal, precision=6)}")
    print(f"resultant: {np.array2string(resultant, precision=4)}")
    print(f"block std rad: {np.array2string(block_std, precision=4)}")
    print(f"accepted: {accepted}")
    print(f"wrote: {args.output}")
    if args.install:
        if not accepted:
            raise RuntimeError("candidate failed quality checks; refusing installation")
        install_calibration(args.output, args.active_path, args.history_dir)
        print(f"installed: {args.active_path}")


if __name__ == "__main__":
    main()
