import numpy as np

from pyEMLearn.core import sqrt, IX, nu, ModeSparse, TransferMatrixSparse,\
                            FieldSparse, KVector, ModeDense, TransferMatrix,\
                            FieldDense
from pyEMLearn.catalog.special import Vacuum
from pyEMLearn.catalog.dielectrics import Air

class Layer:
    _Mode = lambda *args: ModeSparse(*args[1:])
    _Field = lambda *args: FieldSparse(*args[1:])
    _TransferMatrix = TransferMatrixSparse

    def __init__(self,material,thickness):
        self.mat = material
        if self.mat.name is None:
            self.name = self.mat.__class__.__name__
        else:
            self.name = self.mat.name
        self.L = thickness
        self.isFinite = True
        self.isSparse = True

    def __repr__(self):
        return "<Layer: {0}: L = {1}>".format(self.name,self.L)

    def compile(self,gap_layer):
        self.glay = gap_layer
        return self

    def solve(self,kvec):
        self.k = kvec
        self.er = self.mat.er_vec(kvec.wls)
        self.mr = self.mat.mr_vec(kvec.wls)
        self.mode = self._Mode(kvec).solve(self.er,self.mr)
        return self

    def calculate_TM(self):
        # trasfers left gap modes to right gap modes
        self.T = self._TransferMatrix.for_symmetric_layer(
            mode = self.mode,
            mode_gap = self.glay.mode,
            L = self.L
        )
        return self

    def save_field(self,Tsub):
        # Tsub is the transfer matrix from the injection layer to the left gap
        # T_inj_to_inside is the transfer matrix from the injection layer to inside
        self.T_inj_to_inside = Tsub @ self._TransferMatrix.for_interface(
            modeL = self.glay.mode,
            modeR = self.mode
        )
        self.field = self._Field(self.mode)

        return self

    def get_field(self, z, c_inj_p, c_inj_m):
        # z is measured from the left edge of layer
        # c_inj_p and c_inj_m are the mode coefs inside the injection layer
        try:
            c_p,c_m = self.T_inj_to_inside.transfer(c_inj_p, c_inj_m)
            return self.field.get_field_vec(z, c_p, c_m)
        except AttributeError as E:
            raise AttributeError("Internal fields for {} are not saved".format(self.__repr__())) from E

class GapLayer(Layer):
    def __init__(self,material=Vacuum):
        Layer.__init__(self,material,0.0)

    def __repr__(self):
        return "<Gap Layer: {0}>".format(self.name)

    def compile(self,glay):
        raise NotImplementedError
    def save_field(self,Tsub):
        raise NotImplementedError
    def get_field(self, z, c_inj_p, c_inj_m):
        raise NotImplementedError

class HalfInfiniteLayer(Layer):
    def __init__(self,material,injection=True):
        Layer.__init__(self,material,np.inf)
        self.injection = injection
        if self.injection:
            self.direction = "Injection"
        else:
            self.direction = "Transmission"
        self.isFinite = False

    def __repr__(self):
        return "<{0} Layer: {1}>".format(self.direction,self.name)


    def calculate_TM(self):
        if self.injection:
            # transfer injection modes to right gap modes
            self.T = self._TransferMatrix.for_interface(
                modeL = self.mode,
                modeR = self.glay.mode
            )
        else:
            # transfer left gap modes to internal modes for transmission
            self.T = self._TransferMatrix.for_interface(
                modeL = self.glay.mode,
                modeR = self.mode
            )

        return self

    def save_field(self,Tsub):
        # Tsub is the transfer matrix from the injection layer to the left gap
        # T_inj_to_inside is the transfer matrix from the injection layer to inside
        if self.injection:
            self.T_inj_to_inside = TransferMatrixSparse.for_null(2*self.k.kx.size)
        else:
            self.T_inj_to_inside = Tsub @ TransferMatrixSparse.for_interface(
                modeL = self.glay.mode,
                modeR = self.mode
            )
        self.field = FieldSparse(self.mode)

        return self

class LayerPeriodic(Layer):
    _Mode = lambda *args: ModeDense(*args[1:])
    _Field = lambda *args: FieldDense(*args[1:])
    _TransferMatrix = TransferMatrix

    def __init__(self,material,thickness,Lambda_x,Lambda_y):
        Layer.__init__(self,material,thickness)
        self.Lx, self.Ly = Lambda_x, Lambda_y
        self.isSparse = False

    def solve(self,kvec):
        # note kvec must be generated using .from_periodic(...) method

        self.k = kvec
        self.ER_rs = self.mat.er(
            kvec.wls[0],
            kvec.X.flatten(),
            kvec.Y.flatten()
        ).reshape(kvec.X.shape)
        self.MR_rs = self.mat.mr(
            kvec.wls[0],
            kvec.X.flatten(),
            kvec.Y.flatten()
        ).reshape(kvec.X.shape)

        self.ER_fs = np.fft.fftshift(np.fft.fft2(self.ER_rs))/len(kvec.kx) # check
        self.MR_fs = np.fft.fftshift(np.fft.fft2(self.MR_rs))/len(kvec.kx) # check

        S = kvec.X.shape[0]*kvec.X.shape[1]
        self.ER_mat = np.array([[
                    self.ER_fs[tuple(IX(nup,kvec.N,kvec.M)-IX(nupp,kvec.N,kvec.M))]\
                for nupp in range(S)] for nup in range(S)])
        self.MR_mat = np.array([[
                    self.MR_fs[tuple(IX(nup,kvec.N,kvec.M)-IX(nupp,kvec.N,kvec.M))]\
                for nupp in range(S)] for nup in range(S)])

        self.mode = self._Mode(kvec).solve(self.ER_mat,self.MR_mat)

        return self

class System:
    def __init__(self,layer_list=None,
                 injection_material=None,
                 transmission_material=None,
                 gap_material=None):

        if layer_list is None:
            self.layers = []
        else:
            self.layers = layer_list
        self.inj = HalfInfiniteLayer(injection_material or Air,True)
        self.trn = HalfInfiniteLayer(transmission_material or Air,False)
        self.glay = GapLayer(gap_material or Vacuum)

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


    def add_layer(self,layer):
        self.layers += [layer]

    def compile(self):
        self.Lx,self.Ly = [],[]
        self.isPeriodic = False
        for lay in self.layers:
            lay.compile(self.glay)
            if isinstance(lay,LayerPeriodic):
                self.Lx.append(lay.Lx)
                self.Ly.append(lay.Ly)
                self.isPeriodic = True
        assert len(list(set(self.Lx))) < 2 and len(list(set(self.Ly))) < 2,\
            "Cannot compile system with multiple periodic layers of inconsistent dimensions"
        self.inj.compile(self.glay)
        self.trn.compile(self.glay)

        self.left_edges = np.cumsum([0.0]+[l.L for l in self.layers])

        return self

    def solve(self, wls, aoi, save_field = False):
        assert not self.isPeriodic, "System includes periodic layer(s), must use .solve_periodic(...)"
        self.k = KVector.from_AOI(wls, aoi, self.inj.mat)
        self._solve(save_field = save_field)

    def solve_periodic(self, wl, aoi, N, M, save_field = False):
        self.k = KVector.from_periodic(wl,aoi,self.inj.mat,self.Lx[0],self.Ly[0],N,M)
        self._solve(save_field = save_field)

    def _solve(self, save_field):
        self.glay.solve(self.k)

        self.Tcumprod = []

        self.inj.solve(self.k).calculate_TM()
        self.Tcumprod.append(self.inj.T)

        for lay in self.layers:
            lay.solve(self.k).calculate_TM()
            self.Tcumprod.append(self.Tcumprod[-1] @ lay.T)

        self.trn.solve(self.k).calculate_TM()
        self.Tcumprod.append(self.Tcumprod[-1] @ self.trn.T)

        self.S = self.Tcumprod[-1].calculate_SM()

        if save_field:
            self._save_field()

        return self

    def get_incident_mode_coef(self,pol,phase=0):
        # returns c_inj_p the internal mode coef for right waves in injection layer
        # Assumes c_trn_m = 0
        if self.isPeriodic:
            index = len(self.k.kx)//2
        else:
            index = -1
        Etilde_inj_p = self.k.get_incident_field_coef(pol,index=index)*np.exp(1j*phase)
        return self.inj.mode.W_inv @ Etilde_inj_p

    def get_RT(self,pol):
        c_inj_p = self.get_incident_mode_coef(pol)
        # Assumes c_trn_m = 0
        c_inj_m = np.asarray(self.S.S11 @ c_inj_p).flatten()
        c_trn_p = np.asarray(self.S.S21 @ c_inj_p).flatten()

        # print(c_inj_m.shape)

        Etilde_inj_m, Htilde_inj_m = TransferMatrixSparse.for_mode_to_field(
            mode = self.inj.mode
        ).transfer(0*c_inj_p,c_inj_m)
        n = len(Etilde_inj_m)//2
        Ex_inj_m, Ey_inj_m = Etilde_inj_m[:n], Etilde_inj_m[n:]
        Hx_inj_m, Hy_inj_m = Htilde_inj_m[:n], Htilde_inj_m[n:]

        Ez_inj_m = 1j * self.inj.mode._eiky * Hx_inj_m\
                    - 1j * self.inj.mode._eikx * Hy_inj_m

        R = np.abs(Ex_inj_m)**2 + np.abs(Ey_inj_m)**2 + np.abs(Ez_inj_m)**2

        Etilde_trn_p, Htilde_trn_p = TransferMatrixSparse.for_mode_to_field(
            mode = self.trn.mode
        ).transfer(c_trn_p,0*c_trn_p)

        Ex_trn_p, Ey_trn_p = Etilde_trn_p[:n], Etilde_trn_p[n:]
        Hx_trn_p, Hy_trn_p = Htilde_trn_p[:n], Htilde_trn_p[n:]

        Ez_trn_p = 1j * self.trn.mode._eiky * Hx_trn_p\
                    - 1j * self.trn.mode._eikx * Hy_trn_p

        T = np.abs(Ex_trn_p)**2 + np.abs(Ey_trn_p)**2 + np.abs(Ez_trn_p)**2
        kz_trn = sqrt(self.trn.er*self.trn.mr - self.k.kx**2 - self.k.ky**2)
        T = T * (kz_trn/self.trn.mr).real / (self.k.kz_inj/self.inj.mr).real
        # T = T * (self.k.kz_inj/self.inj.mr).real / (kz_trn/self.trn.mr).real

        return R,T

    def get_ellips(self):
        raise NotImplementedError("broken")
        c_inj_p_spol = self.get_incident_mode_coef("s")
        c_inj_m_spol = np.asarray(self.S.S11 @ c_inj_p_spol).flatten()
        c_inj_p_ppol = self.get_incident_mode_coef("p")
        c_inj_m_ppol = np.asarray(self.S.S11 @ c_inj_p_ppol).flatten()
        # rho = (self.Sg.S11 @ np.array([0,1]))[1]/(self.Sg.S11 @ np.array([1,0]))[0]
        # self.PSI = np.arctan(np.abs(rho))*180/np.pi
        # if self.PSI < 0:
        #     self.PSI += 180.0
        # self.DELTA = np.angle(rho,deg=True)
        # if self.DELTA < 0:
        #     self.DELTA += 180.0
        #
        # return self.PSI, self.DELTA

    def _save_field(self):
        self.inj.save_field(None)
        for i,lay in enumerate(self.layers):
            lay.save_field(self.Tcumprod[i])
        self.trn.save_field(self.Tcumprod[-2])

    def get_field(self,z,pol,phase = 0):
        c_inj_p = self.get_incident_mode_coef(pol,phase = phase)
        c_inj_m = np.asarray(self.S.S11 @ c_inj_p).flatten()
        index = np.full(len(z),0,dtype=int)
        for _index,le in enumerate(self.left_edges):
            index[z>le] = _index + 1
        E = np.zeros((index.size,c_inj_m.size//2,3)) + 0j
        mask = index == 0
        if sum(mask) > 0:
            E[mask] = self.inj.get_field(z[mask],c_inj_p,c_inj_m)
        for _index,le in enumerate(self.left_edges):
            mask = (index == _index + 1)
            if sum(mask) > 0:
                E[mask] = self[_index + 1].get_field(z[mask]-le,c_inj_p,c_inj_m)
        return E
