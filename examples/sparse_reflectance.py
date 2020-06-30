import numpy as np
import matplotlib.pyplot as plt

import pyEMLearn.materials as mt
import pyEMLearn.layers as l
from pyEMLearn.catalog.dielectrics import Air, BK7, Pcd2, PMMA
from pyEMLearn.catalog.metals import Ag


s = l.System([
    l.Layer(Ag,0.030),
    l.Layer(Pcd2,0.178),
    l.Layer(Ag,0.030)
]).compile()


wls = np.arange(0.4,0.9,0.005)
aoi = np.arange(0,81,20)
WLS = np.hstack([wls for _ in aoi])
AOI = np.hstack([len(wls)*[_] for _ in aoi])
s.solve(WLS,AOI*np.pi/180)
R,T = s.get_RT("p")
R = R.reshape(len(aoi),len(wls))
T = T.reshape(len(aoi),len(wls))

for RR,TT,a in zip(R,T,aoi):
    plt.plot(wls,RR,label=r"$%.0f\degree$" % a)
    # plt.plot(wls,TT)
    # plt.plot(wls,1-RR-TT)
plt.legend(title="AOI")
plt.ylim(0,1)
plt.xlabel(r"Wavelength $(\mu m)$")
plt.ylabel("Reflectance")
plt.show()




s.solve(0.704,np.pi*0/180,save_field=True)
z = np.arange(-0.3,0.4,0.01)
R,T = s.get_RT("s")

E = s.get_field(z,"s")
P = np.sum(np.abs(E[:,0,:])**2,axis=1)
plt.plot(z,P,color="black",linewidth=2)

plt.ylim(0)
plt.vlines(s.left_edges,*plt.ylim())
plt.text(plt.xlim()[0],plt.ylim()[1],
'''R = %.3f
T = %.3f''' % (R,T),
verticalalignment="top")

plt.xlabel(r"Depth $(\mu m)$")
plt.ylabel(r"Normalized Electric Field Density")

plt.show()
