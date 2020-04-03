#%%
import numpy as np
import matplotlib.pyplot as plt
import pyEMLearn as eml
#%%
Ex = [2.25,2.25+0.16,2.25+0.32]
wls = np.linspace(0.4,0.7,300)
AOI = 0

L1 = 0.05
L2 = 0.072
LA = 0.013

mat1 = eml.materials.CauchyIndex([0.965,0.0221],name="PVA")
mat2 = eml.materials.CauchyIndex([1.026,0.015],name="PMMA")
matA = eml.materials.CauchyIndex([2.7,0],name="PVA")\
        + eml.materials.LorentzIndex(
                A=[0.5,0.4,0.2],
                E=Ex,
                G=[0.1,0.12,0.15]
            )
n1,n2 = mat1.n(0.55)[0], mat2.n(0.55)[0]
Ptot = L1*n1+L2*n2

system = eml.layers.System(
        [
                eml.layers.Layer(
                        eml.catalog.metals.Ag,
                        0.045
                ),
                eml.layers.Layer(
                        mat1,
                        L1
                ),
                eml.layers.Layer(
                        matA,
                        LA
                ),
                eml.layers.Layer(
                        mat2,
                        L2
                ),
                eml.layers.Layer(
                        eml.catalog.metals.Ag,
                        0.045
                )
        ],
        transmission_layer=eml.layers.HalfInfiniteLayer(eml.catalog.dielectrics.BK7,injection=False)
).compile()

X = np.linspace(0.5,1,5)

A = np.zeros((len(X),len(wls)))
for j,x in enumerate(X):
    for i,wl in enumerate(wls):
        system.layers[1].L = (1-x)*Ptot#/n1
        system.layers[3].L = x*Ptot#/n2
        system.solve(wl,AOI*np.pi/180)
        R,T,A[j,i] = system.get_RTA(0)
    print(system)
