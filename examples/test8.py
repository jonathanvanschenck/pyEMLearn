#%%
import numpy as np
import matplotlib.pyplot as plt
import pyEMLearn as eml
#%%

system = eml.layers.System(
        [
                eml.layers.Layer(
                        eml.catalog.metals.Ag,
                        0.045
                ),
                eml.layers.Layer(
                        eml.catalog.dielectrics.Pcd2,
                        0.173
                ),
                eml.layers.Layer(
                        eml.catalog.metals.Ag,
                        0.045
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

plt.figure(figsize=(8,6))
for j,AOI in enumerate(AOIs):
    plt.plot(wls*1e3,A[j],label=r"$%d^\degree$" % int(AOI))#,color=colors[j])
plt.legend(title="AOI")
plt.show()

#%%
#0.52462312, 0.55075377, 0.57286432
wls3 = [0.355,0.45,0.532,0.633]
z = np.linspace(-0.6,0.5,500)
phi = np.linspace(0,2*np.pi*.9999,11)
aveEx = np.zeros((len(wls3),len(z)))
for ii,wl in enumerate(wls3):
    system.solve_fields(wl,AOI = AOIs[0]*np.pi/180)
    Ex = np.zeros((len(z),len(phi)))+0j
    
    for j,_phi in enumerate(phi):
        for i,_z in enumerate(z):
            #print(system.get_internal_field(_z,np.array([1,0])*np.exp(-1j*_phi)))
            Ex[i,j] = system.get_internal_field(_z,np.array([1,0])*np.exp(-1j*_phi))[0]
    #%
    aveEx[ii] = np.mean(Ex.real**2,axis=1)

fig = plt.figure(figsize=(10,6))

#for j,_phi in enumerate(phi):
#    plt.plot(z,Ex[:,j].real**2)
for ii,wl in enumerate(wls3):
    plt.plot(z,aveEx[ii],label=int(wl*1e3))

plt.legend()

#ylim = plt.ylim(-np.max(Ex.real*0)*1.1,np.max(Ex.real**2)*1.1)
ylim = plt.ylim(0,np.max(aveEx)*1.01)

ll = 0
plt.text(ll,ylim[1]*0.9,
         system.inj.name+" ",
         verticalalignment="bottom",
         horizontalalignment="right")
for lay in system.layers:
    plt.axvline(ll,color="black",alpha=0.5)
    plt.text(ll+0.5*lay.L,ylim[1]*0.9,
             lay.name,
             rotation="vertical",
             verticalalignment="bottom",
             horizontalalignment="center")
    ll += lay.L
plt.axvline(ll,color="black",alpha=0.5)
plt.text(ll,ylim[1]*0.9,
         " "+system.trn.name,
         verticalalignment="bottom",
         horizontalalignment="left")
#plt.yscale('log')
#plt.ylim(1e-2)
plt.show()