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

## Storage layout

- `/data1`: raw-voltage Digital RF ring buffer, capped at 10 TB (about 50% of
  the disk).
- `/data2`: analyzed products, including zonal and meridional wind estimates
  under `/data2/products/winds/`.
- `/data3` and `/data4`: available data disks.

## First-cut three-dipole echo selection

`wind_estimates_3ch_2days.py` estimates angle of arrival from the ch1, ch3,
and ch4 dipoles only. It computes normalized cross-spectral coherence on all
three baselines after averaging the auto- and cross-spectra over a 3 × 3
Doppler-range neighborhood. A detection is retained when the mean baseline
coherence is at least 0.80 and its interferometric geometric altitude is
between 70 and 150 km. The AoA search extends to 70 degrees from zenith,
equivalent to 20 degrees elevation above the horizon.

The averaged cross-spectral phases, rather than individual FFT-pixel products,
are passed to the interferometric angle-of-arrival solver. Its distinct
grating-lobe candidates are resolved using the previous accepted
height-resolved mean wind. A lobe is rejected when its observed radial velocity
differs from that wind's predicted radial velocity by more than 40 m/s.

For each accepted AoA, the ch1, ch3, and ch4 complex time series are phase
steered and coherently added. Doppler is then estimated from the centered
one-second complex-voltage segment (20 samples at the 20 Hz coherently
integrated voltage rate). The segment is divided by its RMS amplitude before
the bounded least-squares complex-sinusoid fit; this improves numerical
conditioning without changing fitted frequency or SNR. The hard monostatic
radial-velocity interval is `[-300, 300] m/s`, using
`v_r = -wavelength * f_D / 2`. The fitted sinusoid must also pass the high-SNR
threshold.

The intermediate echo-snippet rows contain timestamp, relative ENU position,
geometric altitude, three direction cosines, slant range, geographic
latitude/longitude/altitude, mean coherence, phase closure, AoA phase residual,
AoA match, fitted Doppler and radial velocity, beamformed SNR, range-gate
index, and mean-wind continuity residual. The timestamp and range-gate index
reference the retained complex RTI snippet in Digital RF metadata.

Wind fits use fitted radial velocities with SNR times squared coherence as
their weight. They are rejected automatically unless there are enough echoes,
at least three occupied azimuth sectors, a well-conditioned inversion, and
acceptable robust velocity residuals. There is no visual or manual acceptance
step; failed bins contain no wind estimate.

Results are written separately under
`/data2/products/winds/3ch_coherent_2day_z70/`. Each window includes zonal and
meridional wind products plus one RTI showing only the automatically selected
pixels used by the wind processor.

## Realtime scheduling

`calc_xc.py` always reduces the newest complete two-second raw-data block
before doing any historical work. It persists a separate backfill cursor and
advances that cursor one block at a time only when the one-minute system load
is below 75% of the available CPU count. Realtime latency is therefore bounded
by at most one in-progress historical block.

The wind timer follows the same policy: it first processes a rolling ten-minute
window ending at the newest reduced two-second sample and maintains a separate
rolling 48-hour realtime product under
`/data2/products/winds/3ch_coherent_realtime_z70/`. Each timer run may process at
most one historical ten-minute block afterward, and only when CPU headroom is
available. Historical plots cannot overwrite the realtime monitor products.
Detections from overlapping realtime windows are de-duplicated by their
two-second timestamp and range-gate index before products and plots are saved.

![rti-1734700851205636](https://github.com/user-attachments/assets/427b2758-fa5a-4433-95e2-4bfb231de57e)
 
![s-1734700848645636](https://github.com/user-attachments/assets/3fb0652c-c6c7-4874-8fba-ff8b9c5b3b52)

![xc-1734700567045636](https://github.com/user-attachments/assets/e0793491-cb5d-42af-95e2-d96cd7c1fb40)
