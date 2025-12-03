import numpy as n
import matplotlib.pyplot as plt
# digitl rf api
import digital_rf as drf
import scipy.signal as ss
import mf_conf as mc
import calc_rti as crti
import os
import time

# read object
d=drf.DigitalRFReader(mc.raw_dir)


dmw = drf.DigitalMetadataWriter(
    mc.xc_dir,
    3600,
    60,
    1000000,
    1,
    "xc",
)
dt=60

dmr = drf.DigitalMetadataReader(mc.xc_dir)
b=d.get_bounds("ch1")

try:
    db=dmr.get_bounds()
    i0=db[1]+dt*1000000
except:
    print("no metadata found")
    i0=b[0]+dt*1000000
print(i0)
# what channels exist
print(d.get_channels())

n_samples=10000000
ipp=10000
n_ipp=int(n_samples/ipp)
offset=8900
dec=10


while True:
    b=d.get_bounds("ch1")
    while (i0 < (b[0]+dt*1000000)):
        print("data no longer exists")
        i0 += dt*1000000
    while ((i0+dt*1000000+1000000) > b[1]):
        print("waiting for more data")
        time.sleep(1)
        b=d.get_bounds("ch1")
    print(i0)
    tvec,rvec,fvec,RTI1,RDI1=crti.rti(d,"ch1",i0,plot=False)

    if False:
        plt.pcolormesh(n.arange(len(tvec)),rvec,n.real(RTI1.T),vmin=-20,vmax=20)
        plt.xlabel("IPP number")
        plt.ylabel("Range (km)")
        cb=plt.colorbar()
        cb.set_label("Real part amplitude")
        plt.show()
    
    for c in ["ch1","ch2","ch3","ch4"]:
        tvec,rvec,fvec,RTI2,RDI2=crti.rti(d,c,i0,plot=False)
#        XC=RTI1.T*n.conj(RTI3.T)
        dB=10.0*n.log10(n.abs(RTI2.T)**2)
        nfloor=n.nanmedian(dB)
        dB=dB-nfloor
        plt.pcolormesh(tvec,rvec,dB,cmap="plasma",vmin=-3,vmax=40)
        cb=plt.colorbar()
        cb.set_label("Power (dB)")
        plt.xlabel("Time (s)")
        plt.ylim([0,400])
        plt.ylabel("Range (km)")
        plt.tight_layout()
        dirname=mc.unix2dirname(i0/1e6)
        dname="/data2/plots/%s/%s"%(c,dirname)
        os.system("mkdir -p %s"%(dname))
        plt.savefig("%s/rti-%06d.png"%(dname,i0/1e6))
        plt.close()
        
    i0+=dt*1000000


    

    
