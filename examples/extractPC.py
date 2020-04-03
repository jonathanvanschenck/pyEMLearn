#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pyEMLearn as eml

#%%

pc = eml.catalog.dielectrics.Pcd2

pc2 = eml.materials.CauchyIndex([1.10946563843109,0.0073377449228674])\
        + eml.materials.LorentzIndex2(
                A = [0.172439520432108,0.0779888680018249,0.0233478966236617,0.0185378206694171,0.0106394246243079],
                E = [1.91598020861984,2.08952480489883,2.25404183731625,2.82841891009979,2.993],
                G = [0.1301859024937,0.0870773365497822,0.0998623211390882,0.0301332046211885,0.0527937552588704]
            )


def get_n(wls,mat):
    return np.array([mat.n(wl)[0] for wl in wls])

def get_k(wls,mat):
    return np.array([mat.n(wl)[1] for wl in wls])

wls = np.arange(0.4,1.3,0.002)


#def diff(p):
#    pc2.matA.coef = np.array(p)
#    return np.sum((get_n(wls,pc)-get_n(wls,pc2))**2+(get_k(wls,pc)-get_k(wls,pc2))**2)
#
#fit = minimize(diff,pc2.matA.coef)
#print(fit)
#pc2.matA.coef = fit.x


fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
ax2 = ax.twinx()


ax.plot(wls,get_n(wls,pc),marker=".",color="red")
ax.plot(wls,get_n(wls,pc2))

ax2.plot(wls,get_k(wls,pc),marker=".",color="green")
ax2.plot(wls,get_k(wls,pc2))