import numpy as n
import matplotlib.pyplot as plt
# digitl rf api
import digital_rf as drf
import mf_conf as mc
# read object
d=drf.DigitalRFReader("/data1/mfraw/")

# what channels exist
print(d.get_channels())

chs=["ch1","ch2","ch3","ch4"]


# bounds for data
b=d.get_bounds("ch1")

for ch in chs:
    nsamp=10000
    offset=900
    z1=d.read_vector_c81d(b[1]-1000000-1300,nsamp,ch)-mc.dc_offset
    print(n.median(z1))
    scale=n.sqrt(n.median(n.abs(z1)**2))
    z1=z1/scale
    plt.figure(figsize=(12,12))
    spec,freq,t,im=plt.specgram(z1,NFFT=256,Fs=1e6,Fc=2.78e6,cmap="plasma",vmin=-80,vmax=20)
    plt.title(ch)
    plt.xlabel("Time (s)")
    plt.ylabel("Frequency (Hz)")
    nfloor=n.median(10.0*n.log10(spec))
    plt.colorbar()
    plt.show()

