#%%
import numpy as np
import matplotlib.pyplot as plt
import pyEMLearn as eml

#from scipy.optimize import minimize
#%%
#wls2 = np.arange(0.6,1.3,0.02)
#n = np.array([eml.catalog.dielectrics.ADTIR.n(wl)[0] for wl in wls2])
#k = np.array([eml.catalog.dielectrics.ADTIR.n(wl)[1] for wl in wls2])
##%%
#
#def cauchy(wls,A,B,C):
#    return A + B/wls**2 + C/wls**4
#
#def diff(p):
#    y = cauchy(wls2,p[0],p[1],p[2])
#    return np.sum((y-n)**2)
#
#fit = minimize(diff,[1.5,0.002,0.002])
##%%
#
#plt.plot(wls2,n,"o")
#plt.plot(wls2,cauchy(wls2,*fit.x))


#%%

s1 = eml.layers.System(
        [
#                eml.layers.Layer(
#                        eml.catalog.metals.Ag,
#                        0.045
#                ),
#                eml.layers.Layer(
#                        eml.catalog.dielectrics.PMMA,
#                        0.00
#                ),
                eml.layers.Layer(
                        eml.catalog.dielectrics.ADTIR,
                        0.05
                ),
                eml.layers.Layer(
                        eml.catalog.dielectrics.PMMA,
                        0.06
                ),
                eml.layers.Layer(
                        eml.catalog.metals.Ag,
                        0.045
                )
        ],
        transmission_layer=eml.layers.HalfInfiniteLayer(eml.catalog.dielectrics.BK7,injection=False)
).compile()
        
        

wls = np.arange(0.5,1.5,0.02)
AOIs = np.arange(20,81,10)

A = np.zeros((len(AOIs),len(wls)))
PSI = np.zeros((len(AOIs),len(wls)))
DEL = np.zeros((len(AOIs),len(wls)))
for j,AOI in enumerate(AOIs):
    for i,wl in enumerate(wls):
        s1.solve(wl,AOI*np.pi/180)
        R,T,A[j,i] = s1.get_RTA(0)
        PSI[j,i],DEL[j,i] = s1.get_ellips()
#%
plt.rcParams['font.size'] = 16
fig = plt.figure(figsize=(8,6))
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
for j,AOI in enumerate(AOIs):
    ax1.plot(wls*1e3,PSI[j],label=r"$%d^\degree$" % int(AOI))
    ax2.plot(wls*1e3,DEL[j],linestyle="--")

plt.rcParams['font.size'] = 12
ax1.legend(title="AOI",loc="lower right")
plt.rcParams['font.size'] = 16
ax1.set_xlabel(r"Wavelength $(nm)$")
ax1.set_ylabel(r"$\Psi$ ",rotation="horizontal")
ax2.set_ylabel(r" $\Delta$",rotation="horizontal")
plt.show()

#%%
plt.rcParams['font.size'] = 16
fig = plt.figure(figsize=(8,6))
ax1 = fig.add_subplot(111)
for j,AOI in enumerate(AOIs):
    ax1.plot(wls*1e3,A[j],label=r"$%d^\degree$" % int(AOI))

plt.rcParams['font.size'] = 12
ax1.legend(title="AOI")
plt.rcParams['font.size'] = 16
ax1.set_xlabel(r"Wavelength $(nm)$")
ax1.set_ylabel(r"A ",rotation="horizontal")
plt.show()