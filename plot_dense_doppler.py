#!/usr/bin/env python3
"""Fit Doppler to every time-range cell in a metadata interval."""

from __future__ import annotations

import argparse
import datetime as dt
import h5py
import os
from pathlib import Path

import digital_rf as drf
import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

import mf_conf as mc
import radar_common as common


CHANNEL_FIELDS = ("rti1", "rti3", "rti4")
READ_CHUNK_S = 60


def parse_utc(value: str) -> dt.datetime:
    parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    return parsed.astimezone(dt.timezone.utc)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("start", type=parse_utc)
    parser.add_argument("end", type=parse_utc)
    parser.add_argument("--metadata-dir", default=mc.xc_dir)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--data-output", type=Path)
    parser.add_argument("--range-min", type=float, default=75.0)
    parser.add_argument("--range-max", type=float, default=125.0)
    parser.add_argument("--display-max-velocity", type=float, default=300.0)
    return parser.parse_args()


def fit_sinusoid_bank(
    times: np.ndarray,
    voltage: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Fit one complex sinusoid to every voltage column.

    The centered one-second segment of every column is normalized to unit RMS.
    A least-squares frequency grid followed by three-point parabolic peak
    interpolation provides a vectorized fit suitable for a dense diagnostic.
    """
    times = np.asarray(times, dtype=np.float64)
    voltage = np.asarray(voltage, dtype=np.complex128)
    if voltage.ndim != 2 or voltage.shape[0] != len(times):
        raise ValueError("voltage must have shape (time, series)")

    sample_interval = float(np.median(np.diff(times)))
    fit_samples = min(
        len(times),
        max(2, int(round(common.DOPPLER_FIT_DURATION_S / sample_interval))),
    )
    fit_start = (len(times) - fit_samples) // 2
    fit_stop = fit_start + fit_samples
    fit_times = times[fit_start:fit_stop] - times[fit_start]
    fit_voltage = voltage[fit_start:fit_stop]

    scale = np.sqrt(np.mean(np.abs(fit_voltage) ** 2, axis=0))
    valid = np.isfinite(scale) & (scale > 0)
    normalized = np.zeros_like(fit_voltage)
    normalized[:, valid] = fit_voltage[:, valid] / scale[valid]

    frequency_grid = np.linspace(
        -common.MAX_FIT_DOPPLER_HZ,
        common.MAX_FIT_DOPPLER_HZ,
        401,
    )
    basis = np.exp(
        -2j * np.pi * frequency_grid[:, None] * fit_times[None, :]
    )
    power = np.abs(basis @ normalized) ** 2
    peak = np.argmax(power, axis=0)
    step = frequency_grid[1] - frequency_grid[0]
    frequency = frequency_grid[peak].copy()

    interior = valid & (peak > 0) & (peak < len(frequency_grid) - 1)
    columns = np.flatnonzero(interior)
    if len(columns):
        center = power[peak[columns], columns]
        lower = power[peak[columns] - 1, columns]
        upper = power[peak[columns] + 1, columns]
        denominator = lower - 2.0 * center + upper
        offset = np.zeros(len(columns), dtype=np.float64)
        curved = np.abs(denominator) > 1e-30
        offset[curved] = 0.5 * (
            lower[curved] - upper[curved]
        ) / denominator[curved]
        frequency[columns] += np.clip(offset, -1.0, 1.0) * step

    model = np.exp(2j * np.pi * fit_times[:, None] * frequency[None, :])
    amplitude = np.mean(np.conj(model) * normalized, axis=0)
    residual = normalized - model * amplitude[None, :]
    residual_power = np.maximum(
        np.mean(np.abs(residual) ** 2, axis=0),
        1e-30,
    )
    snr = np.abs(amplitude) ** 2 / residual_power
    frequency[~valid] = np.nan
    snr[~valid] = np.nan
    return frequency, snr


def fit_sinusoid_fft(
    times: np.ndarray,
    voltage: np.ndarray,
    maximum_frequency_hz: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit a sinusoid to long voltage records using a zero-padded FFT seed."""
    times = np.asarray(times, dtype=np.float64)
    voltage = np.asarray(voltage, dtype=np.complex128)
    if voltage.ndim != 2 or voltage.shape[0] != len(times):
        raise ValueError("voltage must have shape (time, series)")
    if len(times) < 2:
        raise ValueError("at least two voltage samples are required")

    sample_interval = float(np.median(np.diff(times)))
    scale = np.sqrt(np.mean(np.abs(voltage) ** 2, axis=0))
    valid = np.isfinite(scale) & (scale > 0)
    normalized = np.zeros_like(voltage)
    normalized[:, valid] = voltage[:, valid] / scale[valid]

    nfft = max(2048, 2 ** int(np.ceil(np.log2(len(times) * 8))))
    frequency_grid = np.fft.fftshift(
        np.fft.fftfreq(nfft, d=sample_interval)
    )
    spectrum = np.fft.fftshift(
        np.fft.fft(normalized, n=nfft, axis=0),
        axes=0,
    )
    allowed = np.flatnonzero(
        np.abs(frequency_grid) <= maximum_frequency_hz
    )
    allowed_power = np.abs(spectrum[allowed]) ** 2
    peak = np.argmax(allowed_power, axis=0)
    frequency = frequency_grid[allowed[peak]].copy()
    step = frequency_grid[1] - frequency_grid[0]

    interior = valid & (peak > 0) & (peak < len(allowed) - 1)
    columns = np.flatnonzero(interior)
    if len(columns):
        center = allowed_power[peak[columns], columns]
        lower = allowed_power[peak[columns] - 1, columns]
        upper = allowed_power[peak[columns] + 1, columns]
        denominator = lower - 2.0 * center + upper
        offset = np.zeros(len(columns), dtype=np.float64)
        curved = np.abs(denominator) > 1e-30
        offset[curved] = 0.5 * (
            lower[curved] - upper[curved]
        ) / denominator[curved]
        frequency[columns] += np.clip(offset, -1.0, 1.0) * step

    fit_times = times - times[0]
    model = np.exp(2j * np.pi * fit_times[:, None] * frequency[None, :])
    amplitude = np.mean(np.conj(model) * normalized, axis=0)
    residual = normalized - model * amplitude[None, :]
    residual_power = np.maximum(
        np.mean(np.abs(residual) ** 2, axis=0),
        1e-30,
    )
    snr = np.abs(amplitude) ** 2 / residual_power
    frequency[~valid] = np.nan
    snr[~valid] = np.nan
    return frequency, snr


def fit_common_sinusoid_fft(
    times: np.ndarray,
    voltage: np.ndarray,
    maximum_frequency_hz: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Jointly fit one frequency with an independent amplitude per channel.

    ``voltage`` has shape (time, channel, series). Each channel/series is
    normalized to unit RMS, then the four periodogram likelihoods are summed
    before selecting a common frequency for each series.
    """
    times = np.asarray(times, dtype=np.float64)
    voltage = np.asarray(voltage, dtype=np.complex128)
    if voltage.ndim != 3 or voltage.shape[0] != len(times):
        raise ValueError("voltage must have shape (time, channel, series)")
    if len(times) < 2:
        raise ValueError("at least two voltage samples are required")

    sample_interval = float(np.median(np.diff(times)))
    scale = np.sqrt(np.mean(np.abs(voltage) ** 2, axis=0))
    valid_channel = np.isfinite(scale) & (scale > 0)
    valid = np.all(valid_channel, axis=0)
    normalized = np.zeros_like(voltage)
    normalized[:, valid_channel] = (
        voltage[:, valid_channel] / scale[valid_channel]
    )

    nfft = max(2048, 2 ** int(np.ceil(np.log2(len(times) * 8))))
    frequency_grid = np.fft.fftshift(
        np.fft.fftfreq(nfft, d=sample_interval)
    )
    spectrum = np.fft.fftshift(
        np.fft.fft(normalized, n=nfft, axis=0),
        axes=0,
    )
    allowed = np.flatnonzero(
        np.abs(frequency_grid) <= maximum_frequency_hz
    )
    joint_power = np.sum(np.abs(spectrum[allowed]) ** 2, axis=1)
    peak = np.argmax(joint_power, axis=0)
    frequency = frequency_grid[allowed[peak]].copy()
    step = frequency_grid[1] - frequency_grid[0]

    interior = valid & (peak > 0) & (peak < len(allowed) - 1)
    series = np.flatnonzero(interior)
    if len(series):
        center = joint_power[peak[series], series]
        lower = joint_power[peak[series] - 1, series]
        upper = joint_power[peak[series] + 1, series]
        denominator = lower - 2.0 * center + upper
        offset = np.zeros(len(series), dtype=np.float64)
        curved = np.abs(denominator) > 1e-30
        offset[curved] = 0.5 * (
            lower[curved] - upper[curved]
        ) / denominator[curved]
        frequency[series] += np.clip(offset, -1.0, 1.0) * step

    fit_times = times - times[0]
    model = np.exp(2j * np.pi * fit_times[:, None] * frequency[None, :])
    amplitude = np.mean(
        np.conj(model)[:, None, :] * normalized,
        axis=0,
    )
    residual = normalized - model[:, None, :] * amplitude[None, :, :]
    signal_power = np.sum(np.abs(amplitude) ** 2, axis=0)
    residual_power = np.maximum(
        np.sum(np.mean(np.abs(residual) ** 2, axis=0), axis=0),
        1e-30,
    )
    snr = signal_power / residual_power
    fitted_amplitude = amplitude * scale
    frequency[~valid] = np.nan
    snr[~valid] = np.nan
    fitted_amplitude[:, ~valid] = np.nan
    return frequency, snr, fitted_amplitude


def dense_doppler(
    reader: drf.DigitalMetadataReader,
    start_unix: float,
    end_unix: float,
    range_limits: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return time, range, velocity, and selected-channel fit SNR arrays."""
    records = []
    start_us = int(start_unix * 1e6)
    end_us = int(end_unix * 1e6)

    for chunk_start in np.arange(
        start_us,
        end_us,
        int(READ_CHUNK_S * 1e6),
    ):
        chunk_end = min(chunk_start + int(READ_CHUNK_S * 1e6), end_us)
        metadata = reader.read(int(chunk_start), int(chunk_end) - 1)
        for key in sorted(metadata):
            record = metadata[key]
            if not all(
                field in record
                for field in (*CHANNEL_FIELDS, "rvec", "tvec")
            ):
                continue
            record_range = np.asarray(record["rvec"], dtype=np.float64)
            mask = (
                (record_range >= range_limits[0])
                & (record_range <= range_limits[1])
            )
            if not np.any(mask):
                continue
            selected_range = record_range[mask]

            channel_voltage = np.asarray(
                [
                    np.asarray(record[field])[:, mask]
                    for field in CHANNEL_FIELDS
                ],
                dtype=np.complex128,
            )
            sample_count = channel_voltage.shape[1]
            gate_count = channel_voltage.shape[2]
            bank = channel_voltage.transpose(1, 0, 2).reshape(
                sample_count,
                len(CHANNEL_FIELDS) * gate_count,
            )
            frequency, snr = fit_sinusoid_bank(record["tvec"], bank)
            frequency = frequency.reshape(len(CHANNEL_FIELDS), gate_count)
            snr = snr.reshape(len(CHANNEL_FIELDS), gate_count)
            best_channel = np.argmax(
                np.where(np.isfinite(snr), snr, -np.inf),
                axis=0,
            )
            gate = np.arange(gate_count)
            best_frequency = frequency[best_channel, gate]
            best_snr = snr[best_channel, gate]
            records.append(
                (
                    float(key) / 1e6 + 1.0,
                    selected_range,
                    common.DOPPLER_SIGN
                    * common.DOPPLER_TO_MS
                    * best_frequency,
                    best_snr,
                )
            )

    if not records:
        raise RuntimeError("no complete metadata records in requested interval")
    range_km = max((row[1] for row in records), key=len)
    velocity = np.full((len(records), len(range_km)), np.nan)
    fit_snr = np.full_like(velocity, np.nan)
    for row_index, (_, record_range, record_velocity, record_snr) in enumerate(
        records
    ):
        indices = np.searchsorted(range_km, record_range)
        if (
            np.any(indices >= len(range_km))
            or not np.allclose(range_km[indices], record_range)
        ):
            raise ValueError("incompatible range vectors within interval")
        velocity[row_index, indices] = record_velocity
        fit_snr[row_index, indices] = record_snr
    return (
        np.asarray([row[0] for row in records]),
        range_km,
        velocity,
        fit_snr,
    )


def plot_dense_map(
    times_unix: np.ndarray,
    range_km: np.ndarray,
    velocity_ms: np.ndarray,
    start: dt.datetime,
    end: dt.datetime,
    output: Path,
    display_limit_ms: float,
) -> None:
    figure, axis = plt.subplots(
        figsize=(12, 5),
        facecolor="#070b14",
        constrained_layout=True,
    )
    axis.set_facecolor("#050810")
    time_values = [
        dt.datetime.fromtimestamp(value, tz=dt.timezone.utc)
        for value in times_unix
    ]
    image = axis.pcolormesh(
        time_values,
        range_km,
        velocity_ms.T,
        shading="nearest",
        cmap="seismic",
        vmin=-display_limit_ms,
        vmax=display_limit_ms,
    )
    axis.set_xlim(start, end)
    axis.set_ylim(range_km[0], range_km[-1])
    axis.set_xlabel("Time (UTC)", color="#b7c5d9")
    axis.set_ylabel("Round-trip range (km)", color="#b7c5d9")
    axis.set_title(
        "Unfiltered dense Doppler fit · every time–range cell · "
        f"{range_km[0]:.0f}–{range_km[-1]:.0f} km",
        color="#edf4ff",
    )
    axis.tick_params(colors="#8fa1ba")
    axis.xaxis.set_major_formatter(
        mdates.DateFormatter("%H:%M", tz=dt.timezone.utc)
    )
    for spine in axis.spines.values():
        spine.set_color("#314563")
    colorbar = figure.colorbar(image, ax=axis, pad=0.01)
    colorbar.set_label(
        "Monostatic Doppler velocity (m/s)", color="#b7c5d9"
    )
    colorbar.ax.tick_params(colors="#8fa1ba")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".tmp.png")
    figure.savefig(temporary, dpi=180, facecolor=figure.get_facecolor())
    plt.close(figure)
    os.replace(temporary, output)


def main() -> None:
    args = parse_args()
    if args.end <= args.start:
        raise ValueError("end must be after start")
    reader = drf.DigitalMetadataReader(args.metadata_dir)
    values = dense_doppler(
        reader,
        args.start.timestamp(),
        args.end.timestamp(),
        (args.range_min, args.range_max),
    )
    plot_dense_map(
        values[0],
        values[1],
        values[2],
        args.start,
        args.end,
        args.output,
        args.display_max_velocity,
    )
    data_output = args.data_output or args.output.with_suffix(".h5")
    temporary = data_output.with_suffix(data_output.suffix + ".tmp")
    with h5py.File(temporary, "w") as handle:
        handle.create_dataset("time_unix", data=values[0])
        handle.create_dataset("range_km", data=values[1])
        handle.create_dataset(
            "velocity_ms",
            data=values[2],
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
        handle.create_dataset(
            "sinusoid_snr",
            data=values[3],
            compression="gzip",
            compression_opts=1,
            shuffle=True,
        )
    os.replace(temporary, data_output)
    print(f"Plot: {args.output}")
    print(f"Data: {data_output}")
    print(f"Cells: {values[2].size}")


if __name__ == "__main__":
    main()
