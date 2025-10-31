import digital_rf as drf
import matplotlib.pyplot as plt
import numpy as n
import scipy.signal.windows as sw

def lpf(L,om0):
    m=n.arange(L)-L/2.0 + 1e-6
    h=sw.hann(L)*2*n.sin(om0*m)/m
    return(h)

fs=1e6
f1=30e3
om1=(f1/(fs/2))*n.pi
h=lpf(1024,om1)

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


ipp=10000
n_ipp=4000

Z=n.zeros([2,n_ipp,ipp],dtype=n.complex64)
S=n.zeros([2,n_ipp,ipp],dtype=n.float64)

i0=b[1]-n_ipp*ipp-1000000 + 5638
for i in range(n_ipp):
    z0=d.read_vector_c81d(i0+i*ipp,10000,"ch0")
    z1=d.read_vector_c81d(i0+i*ipp,10000,"ch1")

    pha0=n.angle(n.mean(z0[0:80]))
    pha1=n.angle(n.mean(z1[0:80]))

    z0[0:90]=0
    z1[0:90]=0
    z0=n.exp(-1j*pha0)*z0
    z1=n.exp(-1j*pha1)*z1    

    z0f=n.convolve(z0,h,mode="same")
    z1f=n.convolve(z1,h,mode="same")

    Z[0,i,:]=z0f
    Z[1,i,:]=z1f    
    
#    plt.plot(z0.real)
 #   plt.plot(z0.imag)

 #   plt.plot(z0f.real)
#    plt.plot(z0f.imag)    
  #  plt.show()
  #  plt.plot(z1.real)
   # plt.plot(z1.imag)
   # plt.plot(z1f.real)
    #plt.plot(z1f.imag)
    #plt.show()

hw2=sw.hann(n_ipp)
for ri in range(ipp):
    S[0,:,ri]=n.abs(n.fft.fftshift(n.fft.fft(hw2*Z[0,:,ri])))**2.0
    S[1,:,ri]=n.abs(n.fft.fftshift(n.fft.fft(hw2*Z[1,:,ri])))**2.0    

if False:
    plt.pcolormesh(Z.real[0,:,:].T)
    plt.show()

rvec=(n.arange(ipp)-36)*0.15 
fvec=n.fft.fftshift(n.fft.fftfreq(n_ipp,d=ipp/fs))
dB0=10.0*n.log10(S[0,:,:].T)
nfloor=n.nanmedian(dB0)
plt.pcolormesh(fvec,rvec,dB0-nfloor,vmin=0,vmax=50)
plt.ylabel("Range (km)")
plt.xlabel("Doppler shift (Hz)")
plt.colorbar()
plt.show()
dB1=10.0*n.log10(S[1,:,:].T)
nfloor=n.nanmedian(dB1)
plt.pcolormesh(fvec,rvec,dB1-nfloor,vmin=0,vmax=50)
plt.ylabel("Range (km)")
plt.xlabel("Doppler shift (Hz)")

plt.colorbar()
plt.show()
