import numpy as np
import matplotlib.pyplot as plt
import pyEMLearn as eml
#%%
"""
The system is:

Injection Layer |  Layer 1  |  Layer 2   |  Layer 3  | Transmission Layer
 Default (Air)  |   Silver  | Pc (d=2nm) |   Silver  |   BK7 Glass
  L = infty     |  L = 45nm | L = 173nm  |  L = 45nm |   L = infty


We simulate the absorption of this system in the visible range (500-750 nm)
at a series of different angles of incidence (0-80 deg)

"""

system = eml.layers.System(
        [
                eml.layers.Layer( # Layer 1
                        eml.catalog.metals.Ag,
                        0.045 # um
                ),
                eml.layers.Layer( # Layer 2
                        eml.catalog.dielectrics.Pcd2,
                        0.173 # um
                ),
                eml.layers.Layer( # Layer 3
                        eml.catalog.metals.Ag,
                        0.045 # um
                )
        ],
        transmission_layer=eml.layers.HalfInfiniteLayer(eml.catalog.dielectrics.BK7,injection=False)
).compile()

wls = np.linspace(0.5,0.75,100)
AOIs = np.arange(0,81,10)#np.array([0,20,40,60,80])


A = np.zeros((len(AOIs),len(wls)))
PSI = np.zeros((len(AOIs),len(wls)))
DEL = np.zeros((len(AOIs),len(wls)))
for j,AOI in enumerate(AOIs):
    for i,wl in enumerate(wls):
        system.solve(wl,AOI*np.pi/180)
        R,T,A[j,i] = system.get_RTA(90*np.pi/180)
        PSI[j,i],DEL[j,i] = system.get_ellips()

plt.figure()
for i,aoi in enumerate(AOIs):
    plt.plot(wls*1e3,A[i,:],label=r"$%d^{\degree}$" % int(aoi))
plt.legend(title="AOI")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Absorption (1-R-T)")
plt.show()
