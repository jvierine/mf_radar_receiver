import digital_rf as drf
import numpy as n
import matplotlib.pyplot as plt
import mf_conf as mc

dmw = drf.DigitalMetadataWriter(
    "caldir",
    3600,
    3600,
    1000000,
    1,
    "xc",
)

from itertools import combinations

def channel_pairs(ch=["ch1","ch2","ch3","ch4"]):
    chp=[]
    n_ch=len(ch)
    for i in range(n_ch):
        chp.append([i,i])
    idx=n.arange(n_ch)
    chp=chp+list(combinations(idx,2))
    combos=n.array(chp)
    return(combos)

chp=channel_pairs()

dmt = drf.DigitalMetadataReader(mc.xc_dir)
b=dmt.get_bounds()

dt=60

# only use last day
n_t=int(24*3600/dt)
i0=b[1]-24*3600*1000000

noiser0=50
noiser1=110
noisefmin=5
snr_thresh=20



for i in range(n_t):
    print(i)
    dd=dmt.read(i0+i*dt*1000000,i0+(i+1)*dt*1000000)
    xc12=[]
    xc13=[]
    xc14=[]
    
    for k in dd.keys():
        RDI1=dd[k]["rdi1"]
        RDI2=dd[k]["rdi2"]
        RDI3=dd[k]["rdi3"]
        RDI4=dd[k]["rdi4"]
        rdis=[RDI1,RDI2,RDI3,RDI4]
        tvec=dd[k]["tvec"]
        rvec=dd[k]["rvec"]
        fvec=dd[k]["fvec"]
        
        # estimate noise floor
        ri0=n.where(rvec>noiser0)[0][0]
        ri1=n.where(rvec>noiser1)[0][0]
        fi0=n.where(fvec<(-noisefmin))[0][-1]
        fi1=n.where(fvec>noisefmin)[0][0]
        noise_pwr=0.5*(n.mean(n.abs(RDI1[0:fi0,ri0:ri1])**2.0)+n.mean(n.abs(RDI1[fi1:-1,ri0:ri1])**2.0))

        snr=(n.abs(RDI1)**2.0-noise_pwr)/noise_pwr
        rr,ff=n.meshgrid(rvec,fvec)
        idx=n.where( (snr > snr_thresh)&(rvec > 50)&(rvec < 100 )&(n.abs(ff)<1))

        xc12t=(RDI1[idx]*n.conj(RDI2[idx])).flatten()
        xc13t=(RDI1[idx]*n.conj(RDI3[idx])).flatten()
        xc14t=(RDI1[idx]*n.conj(RDI4[idx])).flatten()        
        
        xc12=n.concatenate((xc12,xc12t))
        xc13=n.concatenate((xc13,xc13t))
        xc14=n.concatenate((xc14,xc14t))
    
    odict={"x12":xc12,
           "x13":xc13,
           "x14":xc14}
    dmw.write(i0+i*dt*1000000,odict)



