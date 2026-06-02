import numpy as n
import matplotlib.pyplot as plt
import jcoord
import mf_conf as mc

ngrid=200
u=n.linspace(-1,1,num=ngrid)
v=n.linspace(-1,1,num=ngrid)
uu,vv=n.meshgrid(u,v)
ww=n.sqrt(1-uu**2.0-vv**2.0)

max_zenith_angle=25.0
max_uhor=n.sin(n.pi*max_zenith_angle/180)
idx=n.where( n.sqrt(uu**2 + vv**2.0) < max_uhor)
uu_search=uu[idx]
vv_search=vv[idx]
ww_search=ww[idx]
# this is the search space
uu_search=uu_search.flatten()
vv_search=vv_search.flatten()
ww_search=ww_search.flatten()

# the antennas are in ecef, so we need to convert u,v,w to ecef
ecef=jcoord.enu2ecef(69.58204, 19.22283, 0,uu_search,vv_search,ww_search)
#print(ecef.shape)
#exit(0)


def aoa(xc,posdiff,wavelength=100):
    klen=(2.0*n.pi/wavelength)
    s=n.zeros(len(uu_search),dtype=n.complex64)
    xc2=n.exp(1j*n.angle(xc))
    for i in range(len(xc)):
        s+=xc2[i]*n.exp(1j*klen*(posdiff[i][0]*ecef[0,:]+posdiff[i][1]*ecef[1,:]+posdiff[i][2]*ecef[2,:]))
  #  plt.scatter(uu_search,vv_search,c=n.abs(s)/len(xc),cmap="plasma",vmin=0,vmax=1)
 #   cb=plt.colorbar()
#    cb.set_label("match function")
#    plt.show()
    aoai=n.argmax(n.abs(s))
    return(uu_search[aoai],vv_search[aoai],ww_search[aoai])
