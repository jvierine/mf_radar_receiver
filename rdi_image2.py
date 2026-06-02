import digital_rf as drf
import numpy as n
import matplotlib.pyplot as plt
import mf_conf as mc
import h5py
import os
import jcoord
import image_help as ih

if os.path.exists("phasecal.h5") == False:
    print("no phasecal")
    exit(0)

h=h5py.File("phasecal.h5","r")
phasecal=n.copy(h["phasecal"][()])
h.close()
print(phasecal)
#exit(0)

coord=n.zeros([4,3])
for i in range(4):
    coord[i,:]=jcoord.geodetic2ecef(mc.antenna_coords[i][0],mc.antenna_coords[i][1],mc.antenna_coords[i][2])
    print(coord[i,:])


channels=[0,2,3]

# antenna pos diffs (ignore loop for now)
pos_diffs=[coord[0,:]-coord[2,:],
           coord[0,:]-coord[3,:],
           coord[2,:]-coord[3,:]]

print(coord[0,:]-coord[2,:])
print(pos_diffs[0])
#exit(0)
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
imi=0
for i in range(n_t):
    print(i)
    dd=dmt.read(i0+i*dt*1000000,i0+(i+1)*dt*1000000)
    pts=[]
    for k in dd.keys():
        RDI1=dd[k]["rdi1"]*n.exp(1j*phasecal[0])
        RDI2=dd[k]["rdi2"]*n.exp(1j*phasecal[1])
        RDI3=dd[k]["rdi3"]*n.exp(1j*phasecal[2])
        RDI4=dd[k]["rdi4"]*n.exp(1j*phasecal[3])
#        plt.pcolormesh(10.0*n.log10(n.abs(RDI1)**2.0))
 #       plt.colorbar()
  #      plt.show()

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
        idx=n.where( (snr > snr_thresh)&(rvec > 50)&(rvec < 100 )&(n.abs(ff)<2))
   #     plt.pcolormesh(ff,rr,snr)
    #    plt.colorbar()
     #   plt.show()

        xc12t=(RDI1[idx]*n.conj(RDI2[idx])).flatten()
        xc13t=(RDI1[idx]*n.conj(RDI3[idx])).flatten()
        xc14t=(RDI1[idx]*n.conj(RDI4[idx])).flatten()
        xc34t=(RDI3[idx]*n.conj(RDI4[idx])).flatten()
        rrs=rr[idx].flatten()
        ffs=ff[idx].flatten()
        snrs=snr[idx].flatten()
        
        for mi in range(len(xc12t)):
            u,v,w=ih.aoa([xc13t[mi],xc14t[mi],xc34t[mi]],
                         pos_diffs, wavelength=mc.wavelength)
            print(u,v,w)
            pts.append([rrs[mi]*u,rrs[mi]*v,rrs[mi]*w,ffs[mi],snrs[mi]])
    if len(pts)>0:
        pts=n.array(pts)
        plt.figure(figsize=(10,8))
        plt.subplot(221)
        plt.scatter(pts[:,0],pts[:,1],c=pts[:,3],cmap="brg",vmin=-0.5,vmax=0.5,s=1)
        plt.xlabel("East-West (km)")
        plt.ylabel("North-South (km)")
        plt.colorbar()
        plt.xlim([-60,60])
        plt.ylim([-60,60])

        plt.subplot(222)
        plt.scatter(pts[:,0],pts[:,2],c=pts[:,3],cmap="brg",vmin=-0.5,vmax=0.5,s=1)
        plt.xlabel("East-West (km)")
        plt.ylabel("Up (km)")
        plt.colorbar()
        plt.xlim([-70,70])
        plt.ylim([50,140])

        plt.subplot(223)
        plt.scatter(pts[:,1],pts[:,2],c=pts[:,3],cmap="brg",vmin=-0.5,vmax=0.5,s=1)
        plt.xlabel("North-South (km)")
        plt.ylabel("Up (km)")
        plt.colorbar()
        plt.xlim([-70,70])
        plt.ylim([50,140])

        plt.subplot(224)
        snr[snr<=0]=1e-9
        plt.pcolormesh(fvec,rvec,n.angle(RDI1.T*n.conj(RDI3.T)),cmap="hsv")#,vmin=0,vmax=50)
        plt.ylim([50,200])
        plt.xlim([-2,2])
        plt.xlabel("Doppler (Hz)")
        plt.ylabel("Range (km)")
        plt.colorbar()

            #plt.subplot(222)
            #snr[snr<=0]=1e-9
            #plt.pcolormesh(fvec,rvec,10.0*n.log10(snr.T),cmap="plasma",vmin=0,vmax=30)
            #plt.ylim([50,200])
            #plt.xlim([-2,2])
            #plt.xlabel("Doppler (Hz)")
            #plt.ylabel("Range (km)")
            #plt.colorbar()

        plt.suptitle(mc.unix2datestr((i0+i*dt*1000000)/1e6))
        plt.tight_layout()
        plt.savefig("images/test_img-%06d.png"%(imi))
        plt.close()
        imi+=1
#        plt.show()
                                 
        

