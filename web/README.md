# Ramfjordmoen MF radar monitor

This directory contains the static site published at
`https://juha.no/mf/`.

The receiver runs `plot_monitor_rti.py` to maintain two rolling 48-hour
range-time-intensity plots:

- `latest_rti_48h_full.png`: 0–1500 km one-way range.
- `latest_rti_48h_mesosphere.png`: 0–200 km one-way range.

The plots show quick-look SNR. Five consecutive IPPs and the three dipole
channels are averaged in power, power is reduced to approximately 1.5 km range
bins, and the per-profile noise floor is estimated as the 20th percentile
between 250 and 1400 km. SNR is `(power - noise) / noise`. The mesospheric
color scale is fixed at -3 to 20 dB.

The latest-30-minute 0–200 km RTI instead uses one scalar noise power: the
arithmetic mean power over the entire 30-minute interval and all 30–50 km range
bins. Its 5 IPPs, 10 range samples, and 3 receiver channels provide 150
incoherent power looks. The nominal noise-fluctuation suppression is
`5 log10(150) = 10.9 dB`, so the fixed display scale is -11 to 20 dB. This is
detection gain from reduced variance, not an increase in mean physical SNR.

`publish_monitor.py` collects receiver telemetry and publishes the two plots
and `status.json` to `/var/www/html/mf/` on `juha.no`. The systemd timer in
`systemd/` refreshes the products every five minutes. Static site files are
deployed from this directory and are not overwritten by the receiver.

The realtime wind service also produces a dense, unfiltered fitted-Doppler
RTI for the latest 30 minutes over 50–200 km. Every two-second time–range cell
is fitted independently on each dipole using a centered one-second unit-RMS
complex-voltage segment. The channel with the strongest sinusoid fit is
displayed. The plot uses a fixed monostatic radial-velocity scale of -200 to
+200 m/s; noise-only cells are deliberately retained.
