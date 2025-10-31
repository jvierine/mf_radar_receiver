import digital_rf as drf
import matplotlib.pyplot as plt
import numpy as n
import scipy.signal.windows as sw
import stuffr
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

"""
Signal processing testing for the Ramfjordmoen Medium Frequency Radar
Juha Vierinen, 2024

Block diagram:

RX Dipole 1 ----> |      |
                  | USRP |
RX Dipole 2 ----> |      |

"""
def lpf(L,om0):
    m=n.arange(L)-L/2.0 + 1e-6
    h=sw.hann(L)*2*n.sin(om0*m)/m
    return(h)

# 40 kHz cutoff frequency for the low pass filter
fs=1e6
f1=30e3
om1=(f1/(fs/2))*n.pi
# filter length can't be too long, or there will be ringing of
# strong targets...
h=lpf(256,om1)

show_filter=False
if show_filter:
    plt.plot(h)
    plt.show()
    freqs=n.fft.fftshift(n.fft.fftfreq(len(h),d=1/fs))
    H=n.fft.fftshift(n.abs(n.fft.fft(h))**2.0)
    plt.plot(freqs,10.0*n.log10(H))
    plt.show()
    
d=drf.DigitalRFReader("/data2/mftest_2.78")
print(d.get_channels())
b=d.get_bounds("ch0")
print(b)
i00=b[1]-10e6
z00=d.read_vector_c81d(i00,10000,"ch0")
# estimate BG
#plt.plot(z00.real)
#plt.plot(z00.imag)
#plt.show()
z0bg=n.mean(d.read_vector_c81d(i00,9000000,"ch0"))
z1bg=n.mean(d.read_vector_c81d(i00,9000000,"ch1"))
print(z0bg)
print(z1bg)

ipp=10000
n_ipp=2048
n_avg=1

n_samples = ipp*n_ipp*n_avg
i0=b[0]+1000000+5636
n_t=int(n.floor((b[1]-i0)/n_samples))
print(n_t)
plot_ipp=False

ipp_idx=0
for ti in range(rank,n_t,size):
    Z=n.zeros([2,n_ipp,ipp],dtype=n.complex64)
    S=n.zeros([2,n_ipp,ipp],dtype=n.float64)
    X=n.zeros([1,n_ipp,ipp],dtype=n.complex64)
    cohint=32
    R=n.zeros([3,int(n_ipp/cohint),ipp],dtype=n.complex64)
    tm=n.arange(R.shape[1])*cohint*ipp/fs
    rvec=(n.arange(ipp)-36)*0.15 
    for ai in range(n_avg):
        s0=i0+ti*n_ipp*ipp*n_avg + n_ipp*ipp*ai
        print(ai)
        for i in range(n_ipp):
            z0=d.read_vector_c81d(s0+i*ipp,ipp,"ch0")-z0bg
            z1=d.read_vector_c81d(s0+i*ipp,ipp,"ch1")-z1bg

#            print(z0bg)
 #           print(z1bg)            
  #          z0=z0-z0bg
   #         z1=z1-z1bg
            
            pha0=n.angle(n.mean(z0[0:80]))
            pha1=n.angle(n.mean(z1[0:80]))

            if plot_ipp:
                plt.plot(z0.real)
                plt.plot(z0.imag)        
                plt.show()
#            z0[0:90]=0
 #           z1[0:90]=0
            z0=n.exp(-1j*pha0)*z0
            z1=n.exp(-1j*pha1)*z1    

            z0f=n.convolve(z0,h,mode="same")
            z1f=n.convolve(z1,h,mode="same")

            Z[0,i,:]=z0f
            Z[1,i,:]=z1f
            R[0,int(n.floor(i/cohint)),:]+=z0f
            R[1,int(n.floor(i/cohint)),:]+=z1f
            R[2,int(n.floor(i/cohint)),:]+=z0f*n.conj(z1f)

        hw2=sw.hann(n_ipp)
        dB=10.0*n.log10(n.abs(R[0,:,0:1000])**2.0+n.abs(R[1,:,0:1000])**2.0)
        nfloor=n.median(dB)
        plt.pcolormesh(tm,rvec[0:1000],dB.T-nfloor,vmin=-3,vmax=20,cmap="plasma")
        plt.title(stuffr.unix2datestr(s0/1e6))
        plt.xlabel("Time (s)")
        plt.ylabel("Range (km)")        
        plt.colorbar()
        plt.savefig("rti-%06d.png"%(s0))
        plt.close()

        plt.pcolormesh(tm,rvec[0:1000],n.angle(R[0,:,0:1000].T*n.conj(R[1,:,0:1000].T)),cmap="turbo")
        plt.title(stuffr.unix2datestr(s0/1e6))
        plt.xlabel("Time (s)")
        plt.ylabel("Range (km)")        
        plt.colorbar()
        plt.savefig("arti-%06d.png"%(s0))
        plt.close()

        #        plt.show()
        for ri in range(ipp):
            S0=n.fft.fftshift(n.fft.fft(hw2*Z[0,:,ri]))
            S1=n.fft.fftshift(n.fft.fft(hw2*Z[1,:,ri]))
            S[0,:,ri]+=n.abs(S0)**2.0
            S[1,:,ri]+=n.abs(S1)**2.0
            X[0,:,ri]+=S0*n.conj(S1)


    fvec=n.fft.fftshift(n.fft.fftfreq(n_ipp,d=ipp/fs))
    dB0=10.0*n.log10(S[0,:,:].T)
    nfloor=n.nanmedian(dB0)
    fmax=5
    rmax=200
    dbmax=20
    plt.figure(figsize=(16,10))
    plt.subplot(121)
    plt.pcolormesh(fvec,rvec,dB0-nfloor,vmin=0,vmax=dbmax)
    plt.title("CH0 %s"%(stuffr.unix2datestr(s0/1e6)))
    plt.xlim([-fmax,fmax])
    plt.ylim([0,rmax])    
    plt.ylabel("Range (km)")
    plt.xlabel("Doppler shift (Hz)")
    cb=plt.colorbar()
    cb.set_label("Power (dB)")
#    plt.show()
    dB1=10.0*n.log10(S[1,:,:].T)
    nfloor=n.nanmedian(dB1)
    plt.subplot(122)    
    plt.pcolormesh(fvec,rvec,dB1-nfloor,vmin=0,vmax=dbmax)
    plt.title("CH1")
    plt.xlim([-fmax,fmax])
    plt.ylim([0,rmax])
    plt.ylabel("Range (km)")
    plt.xlabel("Doppler shift (Hz)")
    cb=plt.colorbar()
    cb.set_label("Power (dB)")
    plt.tight_layout()
    plt.savefig("ms-%06d.png"%(s0))
    plt.close()
    
#    plt.show()
    plt.figure(figsize=(16,10))
    plt.subplot(121)

    # remove correlated noise
    S0bg=n.median(S[0,:,:])
    S1bg=n.median(S[1,:,:])
    Xbg=n.median(X[0,:,:])
    X=X-Xbg
    S[0,:,:]=S[0,:,:]-S0bg
    S[1,:,:]=S[1,:,:]-S1bg    
    
    plt.pcolormesh(fvec,rvec,n.angle(X[0,:,:].T),cmap="hsv")
#    plt.xlim([-10,10])
    plt.xlim([-fmax,fmax])
    plt.ylim([0,rmax])    
    
    plt.title("%s"%(stuffr.unix2datestr(s0/1e6)))    
    plt.title("Phase")    
    plt.ylabel("Range (km)")
    plt.xlabel("Doppler shift (Hz)")
    cb=plt.colorbar()
    cb.set_label(r"Phase $\angle \langle z_0(\omega)z_1^*(\omega)\rangle$")
    
    plt.subplot(122)
    
    plt.pcolormesh(fvec,rvec,n.abs(X[0,:,:].T)/(n.sqrt(S[0,:,:].T)*n.sqrt(S[1,:,:].T)),vmin=0,vmax=1)
    plt.xlim([-fmax,fmax])
    plt.ylim([0,rmax])    
    
    plt.xlim([-10,10])    
    plt.title("Coherence")
    plt.ylabel("Range (km)")
    plt.xlabel("Doppler shift (Hz)")
    cb=plt.colorbar()
    cb.set_label("Coherence")
    plt.tight_layout()    
    plt.savefig("mxc-%06d.png"%(s0))
    plt.close()
#    plt.show()
    
