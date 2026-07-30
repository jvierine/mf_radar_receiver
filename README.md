# Software for a new MF radar receiver for Ramfjordmoen and SAURA

Collab between NIPR, UiT, and TGO

## Interferometer phase calibration

`calibrate_interferometer.py` estimates the fixed relative receiver-channel
phases without using the direct transmit path. It processes exact,
non-overlapping ten-second groups using the same joint common-frequency
three-dipole complex-sinusoid fit as the realtime Doppler product. Candidate vertical
echoes must be between 70 and 150 km round-trip range, have fitted Doppler
within ±0.25 Hz, pass 10 dB fitted-sinusoid SNR on all three dipoles, and have
at least 0.80 mean dipole coherence across the five constituent two-second
records.

For each accepted echo, the complex fitted amplitudes provide the channel
phases relative to channel 1. A coherence-weighted circular mean gives the
three-dipole correction. The loop phase is recorded as a diagnostic but is
not calibrated by this method because its phase is not stable over multi-hour
vertical-echo ensembles. For backward compatibility, installation retains
the previous loop entry and marks channel 2 false in
`calibrated_channel_mask`. The HDF5 output contains the correction,
circular concentration and uncertainty, half-hour estimates, sample counts,
selection thresholds, interval, and method. Installation is refused unless
the candidate passes minimum sample-count, circular-concentration, and
block-stability checks. The previous active HDF5 calibration is archived
before an accepted candidate is installed.

Example:

```sh
python calibrate_interferometer.py --hours 24
python calibrate_interferometer.py --hours 24 --install
```

## Synchronized ten-second power and Doppler

`plot_realtime_dense_doppler.py` groups five consecutive two-second metadata
records and fits one common Doppler frequency over the complete ten seconds.
Channels 1, 3, and 4 retain independent complex amplitudes, and their
periodogram likelihoods are summed when selecting the frequency. Every
time-range cell from 0 to 300 km is fitted and retained. No SNR, coherence,
angle-of-arrival, altitude, zenith-angle, or velocity selection is applied.

The HDF5 output contains fitted frequency, monostatic radial velocity,
sinusoid-fit SNR, fitted narrowband power for each dipole, processed broadband
noise power, and the summed narrowband-to-broadband power ratio. The loop is
shown separately for receiver health but excluded from the joint Doppler fit.

The realtime image also contains zonal and meridional winds from an evolving
elliptical full-correlation fit to the three dipoles' power-fading patterns.
Five-minute windows produce one estimate per minute at 70–120 km altitude.
The fitted ground-pattern velocity is divided by two to obtain neutral wind;
correlation and model-residual thresholds blank unreliable estimates without
gating the power or Doppler RTIs.

## Storage layout

- `/data1`: raw-voltage Digital RF ring buffer, capped at 10 TB (about 50% of
  the disk).
- `/data2`: analyzed products and realtime HDF5/plot products.
- `/data3` and `/data4`: available data disks.

## Realtime scheduling

`calc_xc.py` always reduces the newest complete two-second raw-data block
before doing any historical work. It persists a separate backfill cursor and
advances that cursor one block at a time only when the one-minute system load
is below 75% of the available CPU count. Realtime latency is therefore bounded
by at most one in-progress historical block.

The realtime Doppler timer processes the latest 15 minutes into ungated
ten-second fitted-power and Doppler products.

The rolling 48-hour 0–200 km monitor product adds one ungated ten-second
common-frequency Doppler fit per minute below the existing power RTI. Its
incremental HDF5 state is stored on `/data2`.

![rti-1734700851205636](https://github.com/user-attachments/assets/427b2758-fa5a-4433-95e2-4bfb231de57e)
 
![s-1734700848645636](https://github.com/user-attachments/assets/3fb0652c-c6c7-4874-8fba-ff8b9c5b3b52)

![xc-1734700567045636](https://github.com/user-attachments/assets/e0793491-cb5d-42af-95e2-d96cd7c1fb40)
