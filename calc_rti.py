

import numpy as np
import matplotlib.pyplot as plt
import digital_rf as drf
import scipy.constants as sc
import mf_conf as mc


# ── Fixed radar parameters ────────────────────────────────────────────────────
FS_RAW   = 1_000_000   # raw sample rate, Hz
IPP      = 10_000      # samples per IPP at 1 MHz (10 ms)
DEC      = 10          # decimation factor (FIR + downsample)
N_CI     = 5           # IPPs per coherent integration
TX_START = 0           # first sample of TX pulse (relative to IPP start)
TX_LEN   = 100         # samples used to estimate TX phase
GC       = 120         # ground clutter blanking: samples to suppress near TX
OFFSET   = 8900        # sample offset into data block for IPP alignment
TX_CENTER = 54         # TX pulse centre in decimated samples (range offset)
NUM_TAPS  = 50         # FIR filter length

# Range limits
RANGE_MIN_KM = 0.0
RANGE_MAX_KM = 300.0


def range_vector():
    """Decimated range vector (km) for a full IPP."""
    idx  = np.arange(IPP // DEC) * DEC
    rvec = (idx - TX_CENTER) * sc.c / 2 / FS_RAW / 1e3   # km
    return rvec


def range_mask(rvec):
    return (rvec >= RANGE_MIN_KM) & (rvec <= RANGE_MAX_KM)


def rti(d,
        ch="ch1",
        i0=0,
        tx_ch="ch1",
        n_samples=2_000_000,
        fmax=20.0,
        plot=False):
    
    lpf   = mc.fir_lowpass_hann(fc=20e3, fs=FS_RAW, num_taps=NUM_TAPS)
    rvec_full = range_vector()
    rmask     = range_mask(rvec_full)
    rvec      = rvec_full[rmask]
    dec_idx   = np.arange(IPP // DEC) * DEC   # decimated sample indices

    n_ipp    = int(n_samples / IPP)
    n_ci_ipp = n_ipp // N_CI                  # number of output time steps

    n_range = int(rmask.sum())
    RTI = np.zeros((n_ci_ipp, n_range), dtype=np.complex64)
    tvec = np.arange(n_ci_ipp) * N_CI * IPP / FS_RAW

    # ── Step 1: decimate all IPPs into S[n_ipp, n_range_full] ───────────────
    S_full = np.zeros((n_ipp, IPP // DEC), dtype=np.complex64)

    for i in range(n_ipp):
        z    = d.read_vector_c81d(i0 + i * IPP + OFFSET, IPP, ch)    - mc.dc_offset
        z_tx = d.read_vector_c81d(i0 + i * IPP + OFFSET, IPP, tx_ch) - mc.dc_offset

        # TX phase reference from first TX_LEN samples
        tx_phase = np.angle(np.mean(z_tx[:TX_LEN]))

        # Suppress ground clutter region
        z[:GC] = 0.0

        # Phase-correct, filter, decimate
        S_full[i] = np.convolve(
            np.exp(-1j * tx_phase) * z, lpf, mode="same"
        )[dec_idx]

    
    S_ci = S_full[:n_ci_ipp * N_CI].reshape(n_ci_ipp, N_CI, IPP // DEC)
    S_ci = S_ci.sum(axis=1)          # (n_ci_ipp, n_range_full)

    # Keep only the configured round-trip range gates.
    RTI = S_ci[:, rmask].astype(np.complex64)   # (n_ci_ipp, n_range)

    
    prf_ci   = FS_RAW / IPP / N_CI              # effective PRF of CI output
    RDI_full = np.fft.fftshift(
        np.fft.fft(RTI, axis=0), axes=0
    )                                            # (n_ci_ipp, n_range)

    fvec_full = np.fft.fftshift(
        np.fft.fftfreq(n_ci_ipp, d=1.0 / prf_ci)
    )
    fidx = np.where(np.abs(fvec_full) <= fmax)[0]
    fvec = fvec_full[fidx]
    RDI  = RDI_full[fidx, :].astype(np.complex64)   # (n_doppler, n_range)


    return tvec, rvec, fvec, RTI, RDI


# ── Stand-alone test ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    d = drf.DigitalRFReader(mc.raw_dir)
    print("Channels:", d.get_channels())

    b  = d.get_bounds("ch1")
    i0 = b[1] - 2_000_000 - 100_000    # 2 s from end of data

    for ch in ["ch1", "ch2", "ch3", "ch4"]:
        tvec, rvec, fvec, RTI, RDI = rti(d, ch, i0, plot=True)
        print(f"{ch}: RTI {RTI.shape}, RDI {RDI.shape}, "
              f"range {rvec[0]:.1f}–{rvec[-1]:.1f} km, "
              f"Doppler {fvec[0]:.2f}–{fvec[-1]:.2f} Hz")
