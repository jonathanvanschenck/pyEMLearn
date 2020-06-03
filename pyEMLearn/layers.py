"""Docstring for module"""

#%%
import numpy as np

from cmath import sqrt
from numpy.linalg import inv

#from scipy.optimize import minimize

from pyEMLearn.utils import ScatteringMatrix,TransferMatrix
from pyEMLearn.catalog.dielectrics import Air
from pyEMLearn.catalog.special import Vacuum

class Layer:
    """Docstring for Class"""
    def __init__(self,material,thickness):
        """Base class for a layer of material, has finite thickness

        Parameters
        ----------
        material : materials.Material
            An instance of the materials.Material class which describes
            the index of refaction (and permeability) of the layer.
        thickness : float
            The thickness of the layer. The unit of length given here
            must match the unit expectected by materials.Material.n(),
            which is typically um.

        Returns
        -------
        None
        """
        self.mat = material
        if self.mat.name is None:
            self.name = self.mat.__class__.__name__
        else:
            self.name = self.mat.name
        self.L = thickness
        self.isFinite = True


    def __repr__(self):
        return "<Layer: {0}: L = {1}>".format(self.name,self.L)

    def compile(self,gap_layer):
        self.glay = gap_layer

    def solve(self,wl,kx,ky,slow=False):
        """Initalizes layer at wavelenght and kvector (basically AOI)

        Parameters
        ----------
        wl : float
            The (vacuum) wavelength of light passing through the layer.
            The units of length must match that expected by self.mat.n(),
            which is typically um.
        kx : float
            The x componenet of the dimensionless k-vector inside the layer.
        ky : float
            The y componenet of the dimensionless k-vector inside the layer.

        Returns
        -------
        None
        """
        self.wl,self.kx,self.ky = wl,kx,ky
        self.er = self.mat.er(wl)
        self.mr = self.mat.mr(wl)
        ksq = self.er*self.mr

        self.kz = sqrt(ksq-kx**2-ky**2)


        if slow:
            P = np.array([[kx*ky,       ksq-kx**2],
                           [ky**2-ksq,   -kx*ky]])/self.er

            Q = np.array([[kx*ky,     ksq-kx**2],
                          [ky**2-ksq,  -kx*ky]])/self.mr

            Osq = P @ Q

            lamsq, self.W = np.linalg.eig(Osq)
            self.lam = np.vectorize(sqrt)(lamsq)

            self.V = Q @ self.W @ inv(np.diag(self.lam))

            self.W_inv, self.V_inv = inv(self.W), inv(self.V)
        else:
            Q = np.array([[kx*ky,     ksq-kx**2],
                          [ky**2-ksq,  -kx*ky]])/self.mr

            self.lam = np.vectorize(sqrt)(np.array([1j*self.kz,1j*self.kz])**2)#-1j*np.array([self.kz,self.kz])
            self.W = np.eye(2)+0j
            self.W_inv = self.W.copy()
            self.V = Q @ np.diag(1/self.lam)
            self.V_inv = inv(self.V)


    def calculate_SM(self):
        self.S = ScatteringMatrix.for_symmetric_layer(
            W_inv = self.W_inv,
            V_inv = self.V_inv,
            Wg = self.glay.W,
            Vg = self.glay.V,
            lam = self.lam,
            k0 = 2*np.pi/self.wl,
            L = self.L
        )
        return self.S

    def calculate_TM(self,Sg):
        self.T = Sg.get_TM()

    def solve_fields(self):
        self.t = TransferMatrix.for_interface(
            WL = self.glay.W,
            VL = self.glay.V,
            WR_inv = self.W_inv,
            VR_inv = self.V_inv
        )

        self.P = lambda _z: TransferMatrix.for_propegation(
            lam = self.lam,
            k0 = 2*np.pi/self.wl,
            z = _z
        )

        self.C = TransferMatrix.for_conversion(
            W = self.W,
            V = self.V
        )

    def get_internal_field(self,z,c_in_left,c_out_left):
        e = (self.T @ self.t @ self.P(z) @ self.C).transfer(c_in_left,c_out_left)[0]
        e = np.append(e,-(e[0] * self.kx + e[1] * self.ky) / self.kz)
        return e




class GapLayer(Layer):
    """Docstring for Class"""
    def __init__(self,material=Vacuum):
        """Layer class to for gap material

        The gap material is used as to normalized the
        calculations of scattering matricies for regular
        layers in a system. Essentially, the layers are
        symmetrically surrounded by this gap material,
        which has thickness 0.0--which definitionally
        will have no impact on the global optics. But
        we gain simplicity in that all layers see a
        symmetric environment (the gaps) rather than the
        asymmetic environment (the other layers), which
        leads to a much less complicated analytical solution.

        Parameters
        ----------
        material : materials.Material, optional
            An instance of the materials.Material class which describes
            the index of refaction (and permeability) of the gap layer.
            Default is vacuum.

        Returns
        -------
        None
        """
        Layer.__init__(self,material,0.0)

    def __repr__(self):
        return "<Gap Layer: {0}>".format(self.name)

    def solve_fields(self):
        raise NotImplementedError
    def get_internal_field(self,z,c_in_left,c_out_left):
        raise NotImplementedError

class HalfInfiniteLayer(Layer):
    """Docstring for Class"""
    def __init__(self,material,injection=True):
        """Base class for the injection/transmission layers, is assumed
        to semi infinite.

        Parameters
        ----------
        material : materials.Material
            An instance of the materials.Material class which describes
            the index of refaction (and permeability) of the layer.
        injection : bool, optional
            When true, the layer is to be treated as an "injection"
            layer (meaning semi-infinite in the negative direction). If
            false, the layer is treated as a "tranmission" layer (meaning
            semi-infinite in the positive direction). Default is True.

        Returns
        -------
        None
        """
        Layer.__init__(self,material,np.inf)
        self.injection = injection
        if self.injection:
            self.direction = "Injection"
        else:
            self.direction = "Transmission"
        self.isFinite = False

    def __repr__(self):
        return "<{0} Layer: {1}>".format(self.direction,self.name)

    def calculate_SM(self):
        if self.injection:
            self.S = ScatteringMatrix.for_interface(
                WL = self.W,
                VL = self.V,
                WR_inv = self.glay.W_inv,
                VR_inv = self.glay.V_inv
            )
        else:
            self.S = ScatteringMatrix.for_interface(
                WL = self.glay.W,
                VL = self.glay.V,
                WR_inv = self.W_inv,
                VR_inv = self.V_inv
            )
        return self.S

    def calculate_TM(self,Sg=None):
        if self.injection:
            self.T = TransferMatrix(
                np.eye(2) + 0j,
                np.zeros((2,2)) + 0j,
                np.zeros((2,2)) + 0j,
                np.eye(2) + 0j
            )
        else:
            self.T = Sg.get_TM()

    def solve_fields(self):
        self.P = lambda _z: TransferMatrix.for_propegation(
            lam = self.lam,
            k0 = 2*np.pi/self.wl,
            z = _z
        )

        self.C = TransferMatrix.for_conversion(
            W = self.W,
            V = self.V
        )

    def get_internal_field(self,z,c_in_left,c_out_left):
        e = (self.T @ self.P(z) @ self.C).transfer(c_in_left,c_out_left)[0]
        e = np.append(e,-(e[0] * self.kx + e[1] * self.ky) / self.kz)
        return e

class System:
    """Docstring for Class"""
    def __init__(self,layer_list=None,
                 injection_layer=None,
                 transmission_layer=None,
                 gap_layer=None):
        """Base class to describe as solve a system of layers

        The system is composed of a semi-infinite injection layer (from which light is
        incident), a stack of (potentially zero) intermediate layers of finite

        <MORE!>
        """
        if layer_list is None:
            self.layers = []
        else:
            for l in layer_list:
                assert l.isFinite, "Layer Must be finte"
            self.layers = layer_list
        if injection_layer is None:
            self.inj = HalfInfiniteLayer(Air,True)
        else:
            assert not injection_layer.isFinite, "Injection Layer Must be semi-infinte"
            self.inj = injection_layer
        if transmission_layer is None:
            self.trn = HalfInfiniteLayer(Air,False)
        else:
            assert not transmission_layer.isFinite, "Transmission Layer Must be semi-infinte"
            self.trn = transmission_layer
        if gap_layer is None:
            self.glay = GapLayer()
        else:
            self.glay = gap_layer

    def __repr__(self):
        res = "<System:\n"
        res = res + "\t" + self.glay.__repr__() + "\n"
        res = res + "\t" + self.inj.__repr__() + "\n"
        for layer in self.layers:
            res = res + "\t" + layer.__repr__() + "\n"
        res = res + "\t" + self.trn.__repr__() + "\n>"
        return res

    def __getitem__(self,key):
        if type(key) != int:
            raise TypeError
        n = 2 + len(self.layers)
        i = key % n
        if i == 0:
            return self.inj
        elif i == n-1:
            return self.trn
        else:
            return self.layers[i-1]

    def __len__(self):
        return 2 + len(self.layers)
    def __iter__(self):
        self.__i = -1
        return self
    def __next__(self):
        if self.__i < len(self):
            self.__i += 1
            return self[self.__i]
        raise StopIteration


    def add_layer(self,layer,injection=False,transmission=False,gap=False):
        if injection:
            assert not layer.isFinite, "Injection layer must be semi-infinte"
            self.inj = layer
        elif transmission:
            assert not layer.isFinite, "Transmission layer must be semi-infinte"
            self.trn = layer
        elif gap:
            assert layer.L == 0.0, "Gap layer must have zero thickness"
            self.glay = layer
        else:
            assert layer.isFinite, "Layer must be finte"
            self.layers += [layer]

    def compile(self):
        for lay in self.layers:
            lay.compile(self.glay)
        self.inj.compile(self.glay)
        self.trn.compile(self.glay)
        return self

    def solve(self,wl,AOI,save_fields = False,slow=False):
        self.wl, self.AOI = wl, AOI
        n = self.inj.mat.n(wl)
        n = n[0]+n[1]*0j
        kx = n*np.sin(AOI)
        ky = 0.0j

        self.glay.solve(wl,kx,ky,slow=slow)
        self.inj.solve(wl,kx,ky,slow=slow)
        if save_fields:
            self.inj.calculate_TM()
        self.Sg = self.inj.calculate_SM().copy()
        for lay in self.layers:
            lay.solve(wl,kx,ky,slow=slow)
            if save_fields:
                lay.calculate_TM(self.Sg)
            self.Sg = self.Sg @ lay.calculate_SM()
        self.trn.solve(wl,kx,ky,slow=slow)
        self.Sg = self.Sg @ self.trn.calculate_SM()
        if save_fields:
            self.trn.calculate_TM(self.Sg)

        return self

    def get_RTA(self,phase=0):
        pol = np.array([np.sin(phase)*np.cos(self.AOI),np.cos(phase)])
        self.e_inj = self.Sg.S11 @ pol
        E_z = -(self.e_inj[0]*self.inj.kx+self.e_inj[1]*self.inj.ky)/self.inj.kz
        E_ref = np.append(self.e_inj,E_z)
        self.R = np.sum(np.abs(E_ref)**2)
        self.e_trn = self.Sg.S21 @ pol
        E_z = -(self.e_trn[0]*self.trn.kx+self.e_trn[1]*self.trn.ky)/self.trn.kz
        E_trn = np.append(self.e_trn,E_z)
        self.T = np.sum(np.abs(E_trn)**2)\
                     *(self.inj.kz/self.inj.mat.mr(self.wl)).real\
                     /(self.trn.kz/self.trn.mat.mr(self.wl)).real
        self.A = 1-self.T-self.R

        return self.R,self.T,self.A

    def get_ellips(self):
        # Get Ellipsometric Variables
        rho = (self.Sg.S11 @ np.array([0,1]))[1]/(self.Sg.S11 @ np.array([1,0]))[0]
        self.PSI = np.arctan(np.abs(rho))*180/np.pi
        if self.PSI < 0:
            self.PSI += 180.0
        self.DELTA = np.angle(rho,deg=True)
        if self.DELTA < 0:
            self.DELTA += 180.0

        return self.PSI, self.DELTA

    def solve_fields(self,wl,AOI,slow=False):
        self.solve(wl,AOI,save_fields=True,slow=slow)
        self.inj.solve_fields()
        for lay in self.layers:
            lay.solve_fields()
        self.trn.solve_fields()

    def get_internal_field(self,z,pol):
        c_out = self.Sg.S11 @ pol
        zz = 1*z
        if zz < 0:
            return self.inj.get_internal_field(zz,pol,c_out)
        for lay in self.layers:
            if zz < lay.L:
                return lay.get_internal_field(zz,pol,c_out)
            zz -= lay.L
        return self.trn.get_internal_field(zz,pol,c_out)
