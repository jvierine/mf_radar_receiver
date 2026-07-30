# Ramfjordmoen MF radar monitor

This directory contains the static site published at
`https://juha.no/mf/`.

The receiver runs `plot_monitor_rti.py` to maintain two rolling 48-hour
range-time-intensity plots:

- `latest_rti_48h_full.png`: 0–1500 km round-trip range.
- `latest_rti_48h_mesosphere.png`: 0–200 km round-trip range.

The plots show quick-look SNR. Five consecutive IPPs and the three dipole
channels are averaged in power, power is reduced to approximately 1.5 km range
bins, and the per-profile noise floor is estimated as the 20th percentile
between 250 and 1400 km. SNR is `(power - noise) / noise`. The mesospheric
color scale is fixed at -3 to 20 dB.

The upper panel of the latest-15-minute combined product shows the largest
noise-weighted coherent three-dipole matched-filter power on a fixed -20 to
+20 dB scale. Each dipole is referenced to its broadband mean power over the
full 15-minute interval and 30–50 km. Broadband here specifically means
complex-voltage power after transmit-phase correction, the 20 kHz low-pass
filter, 10x decimation, and coherent integration of five 10 ms IPPs. It is not
raw ADC power and is not the sinusoid-fit residual.

`publish_monitor.py` collects receiver telemetry and publishes the two plots
and `status.json` to `/var/www/html/mf/` on `juha.no`. The systemd timer in
`systemd/` refreshes the products every five minutes. Static site files are
deployed from this directory and are not overwritten by the receiver.

The realtime Doppler–AoA service combines synchronized 10-second power and
fitted-Doppler RTIs into one two-panel product for the latest 15 minutes over a
0–300 km display grid. Five consecutive two-second metadata records are
grouped for every output column. At 50–250 km range, the full ten seconds of
phase-calibrated voltage from dipoles 1, 3, and 4 is searched jointly over
Doppler and a 101×101 east/north direction-cosine grid. WGS84 geometry restricts
the search to 50–150 km altitude and at least 20 degrees elevation. Up to
twelve local interferometric ambiguities within 6 dB of the best result are
stored for every time-range cell. The displayed velocity is the strongest
solution and uses a fixed -150 to +150 m/s color scale; noise-only cells are
deliberately retained.

`latest_channel_snr_15m_0_300.png` shows channels 1–4 separately over the
latest 15 minutes. Each column averages the full same 10-second interval used
by the Doppler fit. Each channel is referenced to its own interval-mean power
between 30 and 50 km, while all panels share the same -10 to +20 dB scale.

The old selected-pixel RTI and its realtime processing stage are not used.

`latest_interferometry_ambiguities.png` audits every candidate retained by the
joint Doppler–AoA search. No wind fit or ambiguity selection is currently
applied.
