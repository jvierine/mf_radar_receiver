import numpy as n
import matplotlib.pyplot as plt
# digitl rf api
import digital_rf as drf

# read object
d=drf.DigitalRFReader("/data1/mfraw/")

# what channels exist
print(d.get_channels())

# bounds for data
b=d.get_bounds("ch1")
print(b)

z1=d.read_vector_c81d(b[1]-1000000,40000,"ch1")
z2=d.read_vector_c81d(b[1]-1000000,40000,"ch2")
z3=d.read_vector_c81d(b[1]-1000000,40000,"ch3")
z4=d.read_vector_c81d(b[1]-1000000,40000,"ch4")
plt.figure(figsize=(12,12))
ylim=1e3
plt.subplot(221)
plt.plot(z1.real)
plt.plot(z1.imag)
plt.ylim([-ylim,ylim])
plt.title("CH1")
plt.subplot(222)
plt.plot(z2.real)
plt.plot(z2.imag)
plt.ylim([-ylim,ylim])
plt.title("CH2")
plt.subplot(223)
plt.plot(z3.real)
plt.plot(z3.imag)
plt.ylim([-ylim,ylim])
plt.title("CH3")
plt.subplot(224)
plt.plot(z4.real)
plt.plot(z4.imag)
plt.ylim([-ylim,ylim])
plt.title("CH4")
plt.tight_layout()
plt.show()
print(b)
