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

The upper panel of the latest-15-minute combined product shows
`10 log10(sum_c(|A_c|^2/N_c))` on a fixed -20 to +20 dB scale. `A_c` is the
complex narrowband amplitude at the jointly fitted common Doppler frequency.
`N_c` is the channel's broadband mean power over the full 15-minute interval
and 30–50 km. Broadband here specifically means complex-voltage power after
transmit-phase correction, the 20 kHz low-pass filter, 10x decimation, and
coherent integration of five 10 ms IPPs. It is not raw ADC power and is not
the sinusoid-fit residual.

`publish_monitor.py` collects receiver telemetry and publishes the two plots
and `status.json` to `/var/www/html/mf/` on `juha.no`. The systemd timer in
`systemd/` refreshes the products every five minutes. Static site files are
deployed from this directory and are not overwritten by the receiver.

The realtime wind service combines synchronized 10-second SNR and fitted-Doppler
RTIs into one two-panel product for the latest 15 minutes over 0–300 km.
Five consecutive two-second metadata records are grouped for every output
column. Every time-range cell uses a joint common-frequency fit over the full
10 seconds of unit-RMS complex voltage from channels 1–4. Each channel retains
an independent complex amplitude, and their periodogram likelihoods are summed
before the common Doppler frequency is selected. Doppler uses a fixed
monostatic radial-velocity scale of -200 to +200 m/s; noise-only cells are
deliberately retained.

`latest_channel_snr_15m_0_300.png` shows channels 1–4 separately over the
latest 15 minutes. Each column averages the full same 10-second interval used
by the Doppler fit. Each channel is referenced to its own interval-mean power
between 30 and 50 km, while all panels share the same -10 to +20 dB scale.

The old selected-pixel RTI and its realtime processing stage are not used.
