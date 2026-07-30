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

The upper panel of the latest-30-minute combined product uses one scalar noise
power: the arithmetic mean power over the entire interval and all 30–50 km
range bins. It is sampled every two seconds, matching the dense Doppler panel.
Its 5 IPPs, 10 range samples, and 3 receiver channels provide 150 incoherent
power looks. The nominal noise-fluctuation suppression is
`5 log10(150) = 10.9 dB`, so the fixed display scale is -10 to 20 dB. This is
detection gain from reduced variance, not an increase in mean physical SNR.

`publish_monitor.py` collects receiver telemetry and publishes the two plots
and `status.json` to `/var/www/html/mf/` on `juha.no`. The systemd timer in
`systemd/` refreshes the products every five minutes. Static site files are
deployed from this directory and are not overwritten by the receiver.

The realtime wind service combines the two-second SNR and dense fitted-Doppler
RTIs into one two-panel product for the latest 30 minutes over 0–200 km.
Every time-range cell is fitted independently on each dipole using a centered
one-second unit-RMS complex-voltage segment. The channel with the strongest
sinusoid fit is displayed. Doppler uses a fixed monostatic radial-velocity
scale of -200 to +200 m/s; noise-only cells are deliberately retained.

The old selected-pixel RTI and its realtime processing stage are not used.
