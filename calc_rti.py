import numpy as n
import matplotlib.pyplot as plt
# digitl rf api
import digital_rf as drf
import scipy.signal as ss
import mf_conf as mc
import scipy.constants as sc

def rti(d,
        ch="ch1",
        i0=0,
        tx_ch="ch1",
        gc=120,
        ipp=10000,
        n_samples=10000000,
        offset=8900,
        num_taps=50,
        fmax=10,
        tx_center=54,
        n_window=5, # coherent integration length (how many IPPs)
        plot=False):
    
    dec=10
    lpf=mc.fir_lowpass_hann(fc=20e3, fs=1000000, num_taps=num_taps)
    # how many IPPs to integrate coherently
    tw=n.zeros(n_window)
    tw[:]=1/n_window
    # bounds for data
    idx=n.arange(int(ipp/dec))*dec
    rvec=(idx-tx_center)*sc.c/2/1e6/1e3

    # how many IPPs do we coherently integrate?
    n_ipp=int(n_samples/ipp)
    wf=n.hanning(n_ipp)

    S=n.zeros([n_ipp,int(ipp/dec)],n.complex64)
    tvec=n.arange(n_ipp)*ipp/1e6

    for i in range(n_ipp):
        
        z=d.read_vector_c81d(i0+i*ipp+offset,ipp,ch)-mc.dc_offset
        z_tx=d.read_vector_c81d(i0+i*ipp+offset,ipp,tx_ch)-mc.dc_offset
        if False:
            plt.plot(n.abs(z_tx)**2.0)
            plt.xlabel("Sample index (1 $\mu$s samples)")
            plt.ylabel("Power")
            plt.show()
        
        
        if False:
            plt.plot(z.real)
            plt.plot(z.imag)
            plt.plot(n.abs(z))
            plt.show()
        
        tx_phase=n.angle(n.mean(z_tx[0:100]))
        z[0:gc]=50*z[0:gc]/n.max(n.abs(z[0:gc]))
        S[i,:]=n.convolve(n.exp(-1j*tx_phase)*z,lpf,mode="same")[idx]
        
    nrg=int(ipp/dec)
    FS=n.zeros([n_ipp,int(ipp/dec)],n.complex64)
    for rg in range(nrg):
        FS[:,rg]=n.fft.fftshift(n.fft.fft(wf*S[:,rg]))
        # integrate coherently
        S[:,rg]=n.convolve(tw,S[:,rg],mode="same")

    fvec=n.fft.fftshift(n.fft.fftfreq(n_ipp,d=10e-3))
    fidx=n.where(n.abs(fvec)<fmax)[0]
    
    if plot:
        plt.pcolormesh(tvec,rvec,10.0*n.log10(n.abs(S.T)**2.0),cmap="plasma")
        plt.xlabel("Time (s)")
        plt.ylabel("Range (km)")
        plt.show()

        dB=10.0*n.log10(n.abs(FS.T)**2)
        nfloor=n.median(dB)
        plt.pcolormesh(fvec,rvec,dB,cmap="plasma",vmin=nfloor,vmax=nfloor+50)
        plt.xlim([-5,5])
        plt.colorbar()
        plt.xlabel("Doppler (Hz)")
        plt.ylabel("Range (km)")
        plt.show()
        
    return(tvec,rvec,fvec[fidx],S,FS[fidx,:])

if __name__ == "__main__":
    # read object
    d=drf.DigitalRFReader("/data1/mfraw/")

    # what channels exist
    print(d.get_channels())

    b=d.get_bounds("ch1")
    n_samples=10000000
    i0=b[1]-n_samples-1000000

    for c in  ["ch1","ch2","ch3","ch4"]:
        print(i0)
        tvec,rvec,fvec,RTI,RDI=rti(d,c,i0,plot=True)
