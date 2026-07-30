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

The upper panel of the latest-15-minute combined product uses one scalar noise
power: the arithmetic mean power over the entire interval and all 30–50 km
range bins. Each displayed column uses five consecutive two-second metadata
records, exactly matching the Doppler interval. Each channel contributes 200
power samples after five-IPP coherent integration, and three dipole channels
provide 600 nominal power looks per time/range cell. The fixed display scale
is -10 to 20 dB. This averaging reduces noise fluctuations; it does not
increase the mean physical SNR.

`publish_monitor.py` collects receiver telemetry and publishes the two plots
and `status.json` to `/var/www/html/mf/` on `juha.no`. The systemd timer in
`systemd/` refreshes the products every five minutes. Static site files are
deployed from this directory and are not overwritten by the receiver.

The realtime wind service combines synchronized 10-second SNR and fitted-Doppler
RTIs into one two-panel product for the latest 15 minutes over 0–300 km.
Five consecutive two-second metadata records are grouped for every output
column. Every time-range cell is fitted independently on each dipole using the
full 10 seconds of unit-RMS complex voltage. The channel with the strongest
sinusoid fit is displayed. Doppler uses a fixed monostatic radial-velocity
scale of -200 to +200 m/s; noise-only cells are deliberately retained.

`latest_channel_snr_15m_0_300.png` shows channels 1–4 separately over the
latest 15 minutes. Each column averages the full same 10-second interval used
by the Doppler fit. Each channel is referenced to its own interval-mean power
between 30 and 50 km, while all panels share the same -10 to +20 dB scale.

The old selected-pixel RTI and its realtime processing stage are not used.
