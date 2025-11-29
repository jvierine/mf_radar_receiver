import numpy as n
import matplotlib.pyplot as plt
# digitl rf api
import digital_rf as drf
import scipy.signal as ss

# read object
d=drf.DigitalRFReader("/data1/mfraw/")

# what channels exist
print(d.get_channels())


offset=8900
n_samples=10000000
ipp=10000
# bounds for data
b=d.get_bounds("ch1")

i0=b[1]-n_samples-1000000

n_ipp=int(n_samples/ipp)

dec=10
S=n.zeros([n_ipp,int(ipp/dec)],n.complex64)
ch="ch4"

for i in range(n_ipp):
    z=d.read_vector_c81d(i0+i*ipp+offset,ipp,ch)
    tx_phase=n.angle(n.mean(z[0:100]))
    S[i,:]=ss.decimate(n.exp(-1j*tx_phase)*z,10)
nrg=int(ipp/dec)
FS=n.zeros([n_ipp,int(ipp/dec)],float)
for rg in range(nrg):
    FS[:,rg]=n.fft.fftshift(n.abs(n.fft.fft(S[:,rg]))**2.0)
    
rvec=n.arange(int(ipp/dec))*1.5
tvec=n.arange(n_ipp)*ipp/1e6
plt.pcolormesh(tvec,rvec,10.0*n.log10(n.abs(S.T)**2.0))
plt.xlabel("Time (s)")
plt.ylabel("Range (km)")
plt.show()

fvec=n.fft.fftshift(n.fft.fftfreq(n_ipp,d=10e-3))
dB=10.0*n.log10(FS.T)
nfloor=n.median(dB)
plt.pcolormesh(fvec,rvec,dB,cmap="plasma",vmin=nfloor,vmax=nfloor+50)
plt.colorbar()
plt.xlabel("Doppler (Hz)")
plt.ylabel("Range (km)")
plt.show()
    
