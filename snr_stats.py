import digital_rf as drf
import numpy as n
import matplotlib.pyplot as plt
import mf_conf as mc

dmt = drf.DigitalMetadataReader(mc.xc_dir)
b=dmt.get_bounds()
#print(b)
dd=dmt.read(b[0],b[1])
for k in dd.keys():
    RTI1=n.abs(dd[k]["rti1"])**2.0
    RTI2=n.abs(dd[k]["rti2"])**2.0
    RTI3=n.abs(dd[k]["rti3"])**2.0
    RTI4=n.abs(dd[k]["rti4"])**2.0
    tvec=dd[k]["tvec"]
    rvec=dd[k]["rvec"]
    nfloor1=n.median(RTI1)
    nfloor2=n.median(RTI2)
    SNR1=(RTI1-nfloor1)/nfloor1
    SNR2=(RTI2-nfloor2)/nfloor2
    if False:
        plt.pcolormesh(SNR1.T,vmin=0,vmax=100)
        plt.ylim([0,150])
        plt.colorbar()
        plt.show()
        plt.pcolormesh(SNR2.T,vmin=0,vmax=100)
        plt.ylim([0,150])    
        plt.colorbar()
        plt.show()
    ratio_dB=10.0*n.log10(SNR1.T/SNR2.T)
    ratio_dB[SNR1.T < 10]=n.nan
    plt.figure(figsize=(9,4.8))
    plt.subplot(121)
    ridx=n.where((rvec>50)&(rvec<200))[0]
    print(ratio_dB.shape)
    print(ridx)

    plt.pcolormesh(tvec,rvec[ridx[0]:ridx[-1]],ratio_dB[ridx[0]:ridx[-1],:],vmin=-13,vmax=13,cmap="coolwarm")
    cb=plt.colorbar()
    cb.set_label("Dipole/loop SNR ratio (dB)")
#    plt.ylim([50,200])
    plt.title("SNR ratio (SNR$_1$/SNR$_2$)")
    plt.xlabel("Time (s)")
    plt.ylabel("Range (km)")
    plt.subplot(122)
    plt.hist(ratio_dB[ridx[0]:ridx[-1],:].flatten(),bins=100)
    plt.xlabel("SNR ratio (dB)")
    plt.tight_layout()
    plt.show()
    
#    good_idx=n.where(SNR1.T > 10)
