# Software for a new MF radar receiver for Ramfjordmoen and SAURA

Collab between NIPR, UiT, and TGO

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
between 60 and 150 km.

The averaged cross-spectral phases, rather than individual FFT-pixel products,
are passed to the interferometric angle-of-arrival solver. Wind fits use SNR
times squared coherence as their weight. Results are written separately under
`/data2/products/winds/3ch_coherent_2day/` so they cannot be mixed with older
products that did not have a valid coherence gate.

This is deliberately a first-cut filter. Three long dipole baselines have
grating-lobe ambiguity, and local Doppler-range averaging is only a small
ensemble. The output must be validated against manually inspected intervals
before it is treated as a production wind product.

![rti-1734700851205636](https://github.com/user-attachments/assets/427b2758-fa5a-4433-95e2-4bfb231de57e)
 
![s-1734700848645636](https://github.com/user-attachments/assets/3fb0652c-c6c7-4874-8fba-ff8b9c5b3b52)

![xc-1734700567045636](https://github.com/user-attachments/assets/e0793491-cb5d-42af-95e2-d96cd7c1fb40)
