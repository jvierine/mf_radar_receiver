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
between 250 and 1400 km. The mesospheric color scale is fixed at -3 to 20 dB.

`publish_monitor.py` collects receiver telemetry and publishes the two plots
and `status.json` to `/var/www/html/mf/` on `juha.no`. The systemd timer in
`systemd/` refreshes the products every five minutes. Static site files are
deployed from this directory and are not overwritten by the receiver.
