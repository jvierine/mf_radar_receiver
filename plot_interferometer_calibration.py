import digital_rf as drf
import numpy as n
import matplotlib.pyplot as plt
import mf_conf as mc
import h5py

dmt = drf.DigitalMetadataReader("caldir")
b=dmt.get_bounds()
print(b)

dd=dmt.read(b[0],b[1]-3600*1000000)
xc12=[]
xc13=[]
xc14=[]
phab=n.linspace(-n.pi,n.pi,100)

n_t=0
tk=[]
for k in dd.keys():
    n_t+=1
X12=n.zeros([n_t,99])
X13=n.zeros([n_t,99])
X14=n.zeros([n_t,99])
ti=0
for k in dd.keys():
    print(k)
    tk.append(k)
    xc12=n.concatenate((xc12,dd[k]["x12"]))
    pb,c12=n.histogram(n.angle(dd[k]["x12"]),bins=phab)
    X12[ti,:]=pb#/n.max(pb)
    pb,c12=n.histogram(n.angle(dd[k]["x13"]),bins=phab)
    X13[ti,:]=pb#/n.max(pb)
    pb,c12=n.histogram(n.angle(dd[k]["x14"]),bins=phab)
    X14[ti,:]=pb#/n.max(pb)
    ti+=1
    xc13=n.concatenate((xc13,dd[k]["x13"]))
    
    xc14=n.concatenate((xc14,dd[k]["x14"]))
tk=n.array(tk)
print(X12.shape)
print(len(tk))
plt.subplot(311)
plt.pcolormesh(tk/1e6,phab[0:-1],X12.T)
plt.xlabel("Time (unix)")
plt.ylabel("Phase 1-2 (rad)")
plt.colorbar()
plt.subplot(312)
plt.pcolormesh(tk/1e6,phab[0:-1],X13.T)
plt.ylabel("Phase 1-3 (rad)")
plt.xlabel("Time (unix)")
plt.colorbar()
plt.subplot(313)
plt.pcolormesh(tk/1e6,phab[0:-1],X14.T)
plt.ylabel("Phase 1-4 (rad)")
plt.xlabel("Time (unix)")
plt.colorbar()
plt.tight_layout()
plt.savefig("calpha2.png")
plt.show()

cal_pha=n.zeros(4)
cal_pha[1]=n.angle(n.mean(n.exp(1j*n.angle(xc12))))
cal_pha[2]=n.angle(n.mean(n.exp(1j*n.angle(xc13))))
cal_pha[3]=n.angle(n.mean(n.exp(1j*n.angle(xc14))))
ho=h5py.File("phasecal.h5","w")
ho["phasecal"]=cal_pha
ho.close()
plt.subplot(311)
plt.hist(n.angle(xc12),bins=100)
plt.axvline(n.angle(n.mean(n.exp(1j*n.angle(xc12)))),color="red")
plt.xlabel("Phase 1-2 (rad)")
#plt.axvline(n.angle(n.mean(xc12)))
plt.subplot(312)
plt.hist(n.angle(xc13),bins=100)
plt.axvline(n.angle(n.mean(n.exp(1j*n.angle(xc13)))),color="red")
plt.xlabel("Phase 1-3 (rad)")
#plt.axvline(n.angle(n.mean(xc13)))
plt.subplot(313)
plt.hist(n.angle(xc14),bins=100)
plt.axvline(n.angle(n.mean(n.exp(1j*n.angle(xc14)))),color="red")
plt.xlabel("Phase 1-4 (rad)")
plt.tight_layout()
plt.savefig("calpha.png")
plt.show()
    

