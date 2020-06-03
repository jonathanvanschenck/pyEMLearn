import numpy as np
import matplotlib.pyplot as plt
import pyEMLearn as eml
#%%
L = 0.173
system = eml.layers.System(
        [
                eml.layers.Layer(
                        eml.catalog.metals.Ag,
                        0.045
                ),
                eml.layers.Layer(
                        eml.catalog.special.Vacuum,
                        L
                ),
                eml.layers.Layer(
                        eml.catalog.metals.Ag,
                        0.045
                )
        ],
        transmission_layer=eml.layers.HalfInfiniteLayer(eml.catalog.dielectrics.BK7,injection=False)
).compile()

#%%
lam0 = 1.3*2*L
AOI = 0
system.solve(lam0,AOI*np.pi/180, save_fields=True)
system.solve_fields(lam0,AOI*np.pi/180)
x = np.arange(-0.2,.4,0.01)
y = np.arange(-20,20,0.5)
phi = y*system.inj.kx.real+np.pi*0.5
E = []
for _phi in phi:
    E.append(np.array([system.get_internal_field(_x,[0,np.exp(1j*_phi)]) for _x in x])[:,1])
    plt.plot(x*1e3,E[-1].real**2)
plt.vlines([0,45,45+L*1e3,2*45+L*1e3],*plt.ylim())
plt.show()
E = np.array(E)

X,Y = np.meshgrid(x*1e3,y*1e3)
plt.contourf(X,Y,E.real**2)
plt.vlines([0,45,45+L*1e3,2*45+L*1e3],*plt.ylim())
plt.show()
#%
lam0 = 1.3*2*L
AOI = 10
system.solve(lam0,AOI*np.pi/180, save_fields=True)
system.solve_fields(lam0,AOI*np.pi/180)
x = np.arange(-0.2,.4,0.01)
y = np.arange(-20,20,0.5)
phi = y*system.inj.kx.real+np.pi/2
E = []
for _phi in phi:
    E.append(np.array([system.get_internal_field(_x,[0,np.exp(1j*_phi)]) for _x in x])[:,1])
    plt.plot(x*1e3,E[-1].real**2)
plt.vlines([0,45,45+L*1e3,2*45+L*1e3],*plt.ylim())
plt.show()
E = np.array(E)

X,Y = np.meshgrid(x*1e3,y*1e3)
plt.contourf(X,Y,E.real**2)
plt.vlines([0,45,45+L*1e3,2*45+L*1e3],*plt.ylim())
plt.show()
#%
lam0 = 1.3*2*L
AOI = 20
system.solve(lam0,AOI*np.pi/180, save_fields=True)
system.solve_fields(lam0,AOI*np.pi/180)
x = np.arange(-0.2,.4,0.01)
y = np.arange(-20,20,0.5)
phi = y*system.inj.kx.real+np.pi/2
E = []
for _phi in phi:
    E.append(np.array([system.get_internal_field(_x,[0,np.exp(1j*_phi)]) for _x in x])[:,1])
    plt.plot(x*1e3,E[-1].real**2)
plt.vlines([0,45,45+L*1e3,2*45+L*1e3],*plt.ylim())
plt.show()
E = np.array(E)
#%
X,Y = np.meshgrid(x*1e3,y*1e3)
plt.contourf(X,Y,E.real**2)
plt.vlines([0,45,45+L*1e3,2*45+L*1e3],*plt.ylim())
plt.show()
#%%
lam0 = 1.3*2*L
AOI = 20
system.solve(lam0,AOI*np.pi/180, save_fields=True)
system.solve_fields(lam0,AOI*np.pi/180)
x = np.arange(-0.2,.4,0.01)
y = np.arange(-10,10,0.5)
phi = y*system.inj.kx.real
E = []
for _phi in phi:
    E.append(np.array([system.get_internal_field(_x,[0,np.exp(1j*_phi)]) for _x in x])[:,1])
    plt.plot(x*1e3,E[-1].real**2)
plt.vlines([0,45,45+L*1e3,2*45+L*1e3],*plt.ylim())
plt.show()
E = np.array(E)

X,Y = np.meshgrid(x*1e3,y*1e3)
plt.contourf(X,Y,E.real**2)
plt.vlines([0,45,45+L*1e3,2*45+L*1e3],*plt.ylim())
plt.show()
#%
lam0 = 1.27*2*L
AOI = 20
system.solve(lam0,AOI*np.pi/180, save_fields=True)
system.solve_fields(lam0,AOI*np.pi/180)
x = np.arange(-0.2,.4,0.01)
y = np.arange(-10,10,0.5)
phi = y*system.inj.kx.real
E = []
for _phi in phi:
    E.append(np.array([system.get_internal_field(_x,[0,np.exp(1j*_phi)]) for _x in x])[:,1])
    plt.plot(x*1e3,E[-1].real**2)
plt.vlines([0,45,45+L*1e3,2*45+L*1e3],*plt.ylim())
plt.show()
E = np.array(E)

X,Y = np.meshgrid(x*1e3,y*1e3)
plt.contourf(X,Y,E.real**2)
plt.vlines([0,45,45+L*1e3,2*45+L*1e3],*plt.ylim())
plt.show()
#%
lam0 = 1.25*2*L
AOI = 20
system.solve(lam0,AOI*np.pi/180, save_fields=True)
system.solve_fields(lam0,AOI*np.pi/180)
x = np.arange(-0.2,.4,0.01)
y = np.arange(-10,10,0.5)
phi = y*system.inj.kx.real
E = []
for _phi in phi:
    E.append(np.array([system.get_internal_field(_x,[0,np.exp(1j*_phi)]) for _x in x])[:,1])
    plt.plot(x*1e3,E[-1].real**2)
plt.vlines([0,45,45+L*1e3,2*45+L*1e3],*plt.ylim())
plt.show()
E = np.array(E)
#%
X,Y = np.meshgrid(x*1e3,y*1e3)
plt.contourf(X,Y,E.real**2)
plt.vlines([0,45,45+L*1e3,2*45+L*1e3],*plt.ylim())
plt.show()