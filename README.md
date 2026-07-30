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

## Joint Doppler and angle-of-arrival fit

`aoa_doppler.py` searches each complete ten-second block jointly over Doppler
and a regular east/north direction-cosine grid using only the calibrated
dipoles on channels 1, 3, and 4. At each range from 50 to 250 km, `jcoord`
converts every trial direction to WGS84 geodetic altitude. Only directions at
50–150 km altitude and at least 20 degrees elevation are searched.

The largest coherent three-dipole matched-filter value supplies the displayed
signal power and radial velocity. Distinct local direction-grid maxima within
6 dB of the best solution are retained rather than unwrapped. Doppler maxima
are tested in descending incoherent power. That power is an upper bound on
every coherent beam, so the search stops only when no remaining Doppler can
fall within the retained 6 dB interval. The realtime HDF5 product stores up to
twelve ambiguities per time-range cell, including direction cosines, WGS84
position and altitude, Doppler, radial velocity, coherent power, relative
power, and phase match.

## Storage layout

- `/data1`: raw-voltage Digital RF ring buffer, capped at 10 TB (about 50% of
  the disk).
- `/data2`: analyzed products and realtime HDF5/plot products.
- `/data3` and `/data4`: available data disks.

## Interferometry diagnostics

`plot_interferometry_debug.py` reads the ambiguity table produced by the joint
ten-second matched filter. It writes an ambiguity audit and a deliberately
labelled provisional zonal/meridional fit. The diagnostic alternates between
selecting one Doppler–AoA candidate per time-range cell and fitting zonal and
meridional wind on 10 km altitude knots with second-difference smoothness.
Only candidates with at least 10 dB coherent-power ratio, 0.80 phase match,
70–150 km altitude, and radial velocity within ±150 m/s enter this diagnostic.
A 40 m/s residual gate is applied only when displaying the provisional fit.

The selected candidates, predicted radial velocities, residuals, acceptance
mask, and fitted profiles are stored in
`latest_interferometry_debug.h5`. This is a debugging product, not yet a
validated wind retrieval.

## Realtime scheduling

`calc_xc.py` always reduces the newest complete two-second raw-data block
before doing any historical work. It persists a separate backfill cursor and
advances that cursor one block at a time only when the one-minute system load
is below 75% of the available CPU count. Realtime latency is therefore bounded
by at most one in-progress historical block.

The realtime interferometry timer processes the latest 15 minutes into
ten-second joint Doppler–AoA fits, writes all retained ambiguities to HDF5, and
then refreshes the two diagnostic plots. It does not run the removed
selected-pixel wind processor.

![rti-1734700851205636](https://github.com/user-attachments/assets/427b2758-fa5a-4433-95e2-4bfb231de57e)
 
![s-1734700848645636](https://github.com/user-attachments/assets/3fb0652c-c6c7-4874-8fba-ff8b9c5b3b52)

![xc-1734700567045636](https://github.com/user-attachments/assets/e0793491-cb5d-42af-95e2-d96cd7c1fb40)
