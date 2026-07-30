# Ramfjordmoen MF radar monitor

This directory contains the static site published at
`https://juha.no/mf/`.

The receiver runs `plot_monitor_rti.py` to maintain two rolling 48-hour
range-time-intensity plots:

- `latest_rti_48h_full.png`: 0–1500 km round-trip range.
- `latest_rti_48h_mesosphere.png`: synchronized noise-referenced power and
  fitted Doppler panels over 0–200 km round-trip range.

The plots show quick-look SNR. Five consecutive IPPs and the three dipole
channels are averaged in power, power is reduced to approximately 1.5 km range
bins, and the per-profile noise floor is estimated as the 20th percentile
between 250 and 1400 km. SNR is `(power - noise) / noise`. The mesospheric
color scale is fixed at -3 to 20 dB.

The upper panel of the latest-15-minute combined product shows the sum of the
three independently fitted dipole powers on a fixed -20 to +20 dB scale. Each
dipole is referenced to its broadband mean power over the full 15-minute
interval and 30–50 km. Broadband here specifically means
complex-voltage power after transmit-phase correction, the 20 kHz low-pass
filter, 10x decimation, and coherent integration of five 10 ms IPPs. It is not
raw ADC power and is not the sinusoid-fit residual.

`publish_monitor.py` collects receiver telemetry and publishes the plots
and `status.json` to `/var/www/html/mf/` on `juha.no`. The systemd timer in
`systemd/` refreshes the products every five minutes. Static site files are
deployed from this directory and are not overwritten by the receiver.

The realtime Doppler service combines synchronized 10-second fitted-power and
Doppler RTIs with two fading-wind panels for the latest 15 minutes. The RTIs
cover 0–300 km. Five consecutive two-second metadata records are grouped for
every output column. Every time-range cell is fitted. The three dipoles share
one Doppler frequency but retain independent complex amplitudes. No SNR,
coherence, AoA, altitude, zenith-angle, or velocity gate is applied.

The same image includes zonal and meridional fading-wind panels at 70–120 km.
Each minute uses a five-minute window of power fluctuations from the three
dipoles. An evolving elliptical correlation model supplies the ground-pattern
velocity, which is divided by two for neutral wind. Wind cells are retained
only when median cross-correlation is at least 0.15 and normalized model RMSE
is at most 0.35; these automated quality gates do not alter either RTI.

`latest_fading_correlation_diagnostics_15m.png` shows the unmasked
cross-correlation peak delay and peak coefficient for each of the three
dipole baselines. This makes the raw information supporting—or contradicting—
each fading-wind estimate directly visible.

`latest_channel_snr_15m_0_300.png` shows channels 1–4 separately over the
latest 15 minutes. Each column averages the full same 10-second interval used
by the Doppler fit. Each channel is referenced to its own interval-mean power
between 30 and 50 km, while all panels share the same -10 to +20 dB scale.

The old selected-pixel RTI and its realtime processing stage are not used.

The 48-hour mesospheric product retains one ten-second common-frequency
three-dipole Doppler fit per minute in `doppler_48h_state.h5`. Every available
0–200 km cell is shown without SNR, coherence, altitude, or echo-classification
selection. The Doppler color scale is fixed at -150 to +150 m/s.
