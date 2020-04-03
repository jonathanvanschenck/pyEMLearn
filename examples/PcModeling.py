#%%
import numpy as np
import matplotlib.pyplot as plt
import pyEMLearn as eml

#%%

L = 0.26#0.165
Lag = 0.040

BK7 = eml.catalog.dielectrics.BK7
Ag = eml.catalog.metals.Ag
Pc = eml.catalog.dielectrics.Pcd2



film = eml.layers.System(
        [
                eml.layers.Layer(Pc,L)
        ],
        injection_layer=eml.layers.HalfInfiniteLayer(BK7,injection=True)
).compile()

bottom = eml.layers.System(
        [
                eml.layers.Layer(Ag,Lag),
                eml.layers.Layer(Pc,L)
        ],
        injection_layer=eml.layers.HalfInfiniteLayer(BK7,injection=True)
).compile()

top = eml.layers.System(
        [
                eml.layers.Layer(Pc,L),
                eml.layers.Layer(Ag,Lag)
        ],
        injection_layer=eml.layers.HalfInfiniteLayer(BK7,injection=True)
).compile()

cavity = eml.layers.System(
        [
                eml.layers.Layer(Ag,Lag),
                eml.layers.Layer(Pc,L),
                eml.layers.Layer(Ag,Lag)
        ],
        injection_layer=eml.layers.HalfInfiniteLayer(BK7,injection=True)
).compile()


sys = [film,bottom,top,cavity][-2:]
sysN = ["film","bottom","top","cavity"][-2:]
#%

wls = np.arange(0.4,1.0,0.001)
AOI = 0

R = np.zeros((len(sys),len(wls)))
T = np.zeros((len(sys),len(wls)))
A = np.zeros((len(sys),len(wls)))
for j,s in enumerate(sys):
    for i,wl in enumerate(wls):
        s.solve(wl,AOI*np.pi/180)
        R[j,i],T[j,i],A[j,i] = s.get_RTA(0)
        
#%
plt.figure(figsize=(8,6))
for j,s in enumerate(sys):
    plt.plot(wls,R[j],label=sysN[j])
plt.axvline(.633,color="black",linestyle="--")
plt.legend()
plt.ylim(0)
plt.show()

#%

wl = 0.633
z = np.linspace(-0.2,0.5,200)
phi = np.linspace(0,2*np.pi*.9999,11)
aveEx = np.zeros((len(sys),len(z)))

for ii,s in enumerate(sys):
    s.solve_fields(wl,AOI = AOI*np.pi/180)
    Ex = np.zeros((len(z),len(phi)))+0j
    
    for j,_phi in enumerate(phi):
        for i,_z in enumerate(z):
            #print(system.get_internal_field(_z,np.array([1,0])*np.exp(-1j*_phi)))
            Ex[i,j] = s.get_internal_field(_z,np.array([1,0])*np.exp(-1j*_phi))[0]
    #%
    aveEx[ii] = np.mean(Ex.real**2,axis=1)
    
    #%
    
colors = ["black","red","green","blue"]
fig = plt.figure(figsize=(10,6))

for ii,s in enumerate(sys):
    
    plt.plot(z,aveEx[ii],label=sysN[ii],color=colors[ii])
    ll = 0
    for lay in s.layers:
        plt.axvline(ll,alpha=0.5,color=colors[ii])
        ll += lay.L
    plt.axvline(ll,color=colors[ii],alpha=0.5)
    
plt.legend()
ylim = plt.ylim(0,np.max(aveEx)*1.01)
plt.show()