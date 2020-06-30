import numpy as np
from scipy import sparse as sp
from scipy.sparse import linalg as sla

from cmath import sqrt as _sqrt
sqrt = np.vectorize(_sqrt)

def cross(v1,v2):
    return np.array([
        v1[1]*v2[2]-v1[2]*v2[1],
        v1[2]*v2[0]-v1[0]*v2[2],
        v1[0]*v2[1]-v1[1]*v2[0]
    ])


def norm_cross(vec1,vec2):
    result = cross(vec1,vec2)*(1+0j)
    return result/np.linalg.norm(result)
norm_cross_vec = np.vectorize(norm_cross,signature='(n),(n)->(n)')

def get_s_pol(vec):
    result = np.cross([0j,0j,1+0j],vec)
    l = np.linalg.norm(result)
    if l == 0.0:
        return np.array([1+0j,0j,0j])
    return result/l
get_s_pol_vec = np.vectorize(get_s_pol,signature='(n)->(n)')


def IX(nu,N,M):
    return np.array([(nu//(2*M+1)), nu %(2*M+1)])
def nu(ni,mi,N,M):
    return (ni%(2*N+1))*(2*M+1) + (mi%(2*M+1))

def linspace(lb,rb,num):
    if num == 1:
        return np.array([(rb+lb)/2])
    return np.linspace(lb,rb,num)

class KVector:
    def __init__(self,kx,ky,ksq,wls):
        self.kx,self.ky,self.ksq,self.wls = kx,ky,ksq,wls
        self.kz_inj = sqrt(self.ksq - self.kx**2 - self.ky**2)
        k0 = np.vstack([self.kx, self.ky, self.kz_inj]).T
        s_pol = get_s_pol_vec(k0)
        p_pol = norm_cross_vec(k0, s_pol)
        self.pol_vec = {"s":s_pol, "p":p_pol}

    def __repr__(self):
        return "<KVector Object>"

    def get_kz(self,material):
        raise NotImplementedError("Do I need this?")

    def get_incident_field_coef(self,pol,index=-1):
        # returns the Etilde = [Ex,Ey] vector for the incident field
        # pol is a string for either "s" or "p" polarized light
        # index specifies which kvec to populate, if -1 it populates all kvects
        ExEy = np.hstack(self.pol_vec[pol].T[:2])
        if index == -1:
            return ExEy
        ExEy_sub = np.zeros(ExEy.shape) + 0j
        ExEy_sub[index] = ExEy[index]
        ExEy_sub[2*index] = ExEy[2*index]
        return ExEy_sub

    @classmethod
    def from_AOI(cls,wls,aoi,inj_mat):
        WLS = np.atleast_1d(wls)
        AOI = np.atleast_1d(aoi)
        ER = np.vectorize(inj_mat.er)(WLS)
        MR = np.vectorize(inj_mat.mr)(WLS)
        ksq = ER*MR
        ktot = sqrt(ksq)
        ky = np.sin(AOI)*ktot
        new_cls = cls(0*ky,ky,ksq,wls)
        new_cls.aoi = AOI
        return new_cls

    @classmethod
    def from_periodic(cls,wl,aoi,inj_mat,Lambda_x,Lambda_y,N,M):
        S = (2*N+1)*(2*M+1)
        tx, ty = wl/Lambda_x,wl/Lambda_y#2*np.pi/Lambda_x, 2*np.pi/Lambda_y
        ksq = inj_mat.er(wl)*inj_mat.mr(wl)
        ky = np.sin(aoi)*_sqrt(ksq)
        kz = np.cos(aoi)*_sqrt(ksq)

        Kx = 0 + (np.arange(S) // (2*M+1) - N)*tx
        Ky = ky + (np.arange(S) % (2*M+1) - M)*ty
        Ksq = Kx**2 + Ky**2 + kz**2

        new_cls = cls(
            Kx,
            Ky,
            Ksq,
            np.full(S,wl)
        )

        new_cls.N, new_cls.M = N,M
        new_cls.Lambda_x, new_cls.Lambda_y = Lambda_x, Lambda_y
        # X is a 2N+1x2M+1 dim vector of x-vals: [[x0,x0,...],[x1,x1,...],...]
        # Y is a 2N+1x2M+1 dim matrix of y-vals: [[y0,y1,...],[y0,y1,...],...]
        new_cls.Y,new_cls.X = np.meshgrid(linspace(0,Lambda_y,(2*M+1)),
                                          linspace(0,Lambda_x,(2*N+1)))

        return new_cls

class Mode:
    def __init__(self,kvec):
        self.k = kvec

    def solve(self,er,mr):
        raise NotImplementedError("Virtual Method")

class ModeDense(Mode):
    def solve(self,ER_mat,MR_mat):
        ER_mat_inv = np.linalg.inv(ER_mat)
        MR_mat_inv = np.linalg.inv(MR_mat)

        Kx_mat = np.diag(self.k.kx)
        Ky_mat = np.diag(self.k.ky)

        eix = ER_mat_inv @ Kx_mat
        eiy = ER_mat_inv @ Ky_mat
        mix = MR_mat_inv @ Kx_mat
        miy = MR_mat_inv @ Ky_mat

        # Need for Ez calc in field
        self._eikx = eix
        self._eiky = eiy

        eixx = Kx_mat @ eix
        eiyx = Ky_mat @ eix
        eixy = Kx_mat @ eiy
        eiyy = Ky_mat @ eiy

        mixx = Kx_mat @ mix
        miyx = Ky_mat @ mix
        mixy = Kx_mat @ miy
        miyy = Ky_mat @ miy

        self.P = np.bmat([[eixy,MR_mat-eixx],[eiyy-MR_mat,-eiyx]])
        self.Q = np.bmat([[mixy,ER_mat-mixx],[miyy-ER_mat,-miyx]])
        Osq = self.P @ self.Q
        lamsq, self.W = np.linalg.eig(Osq)
        self.lam = sqrt(lamsq)
        self.X = lambda zp: np.diag(np.exp(self.lam * zp))
        self.W_inv = self.W.T.conj()
        self.V = self.Q @ self.W @ np.diag(1/self.lam)
        self.V_inv = np.diag(1/self.lam) @ self.W_inv @ self.P

        return self


class ModeSparse(Mode):
    def solve(self, er, mr):
        ks = er * mr
        kxs, kys, kxy = self.k.kx**2, self.k.ky**2 , self.k.kx*self.k.ky
        kz = sqrt(ks - kxs - kys)

        # Need for Ez calc in field
        self._eikx = self.k.kx/er
        self._eiky = self.k.ky/er

        self.P = sp.bmat([
            [         sp.diags(kxy/er,format="csc"), sp.diags(ks/er - kxs/er,format="csc") ],
            [ sp.diags(kys/er - ks/er,format="csc"),        -sp.diags(kxy/er,format="csc") ]
        ],format="csc")
        self.Q = sp.bmat([
            [         sp.diags(kxy/mr,format="csc"), sp.diags(ks/mr - kxs/mr,format="csc") ],
            [ sp.diags(kys/mr - ks/mr,format="csc"),        -sp.diags(kxy/mr,format="csc") ]
        ],format="csc")

        _lam = sqrt((1j*kz)**2)
        self.lam = np.append(_lam,_lam)
        # self.X = lambda zp: np.diag(np.exp(self.lam * zp))

        self.W = sp.eye(len(self.lam),format="csc")
        self.W_inv = self.W.copy()
        self.V = self.Q @ sp.diags(1/self.lam,format="csc")
        self.V_inv = sp.diags(1/self.lam,format="csc") @ self.P

        return self

class Field:
    def __init__(self,mode):
        self.mode = mode
        N = mode.k.kx.size
        # signature = "(n),({0}),({0})->({1},3)".format(2*N,N)
        signature = "(),({0}),({0})->({1},3)".format(2*N,N)
        self.get_field_vec = np.vectorize(self.get_field,signature=signature)

    def get_field(self,z,c_p,c_m):
        raise NotImplementedError("Virtual Method")

class FieldDense(Field):
    def __init__(self,*arg,**kwarg):
        Field.__init__(self,*arg,**kwarg)
        self.T_mode_to_field = TransferMatrix.for_mode_to_field(self.mode)

    def get_field(self,z,c_p,c_m):
        # c_p and c_m are internal mode
        zp = 2*np.pi*z/self.mode.k.wls
        T = TransferMatrix.for_propegation(self.mode.lam,zp) @ self.T_mode_to_field
        Etilde,Htilde = T.transfer(c_p,c_m)
        n = len(Etilde)//2
        Ex,Ey = Etilde[:n], Etilde[n:]
        Hx,Hy = Htilde[:n], Htilde[n:]
        Ez = 1j*self.mode._eiky @ Hx - 1j*self.mode._eikx @ Hy
        return np.vstack([Ex,Ey,Ez]).T

class FieldSparse(Field):
    def __init__(self,*arg,**kwarg):
        Field.__init__(self,*arg,**kwarg)
        self.T_mode_to_field = TransferMatrixSparse.for_mode_to_field(self.mode)

    def get_field(self,z,c_p,c_m):
        # c_p and c_m are internal modes
        zp = 2*np.pi*z/self.mode.k.wls
        T = TransferMatrixSparse.for_propegation(self.mode.lam,zp) @ self.T_mode_to_field
        Etilde,Htilde = T.transfer(c_p,c_m)
        n = len(Etilde)//2
        Ex,Ey = Etilde[:n], Etilde[n:]
        Hx,Hy = Htilde[:n], Htilde[n:]
        Ez = 1j * self.mode._eiky * Hx - 1j * self.mode._eikx * Hy
        return np.vstack([Ex,Ey,Ez]).T


class TransferMatrix:
    _inv = lambda *arg: np.linalg.inv(arg[-1])
    _sqzeros = lambda *arg: np.zeros((arg[-1],arg[-1]))
    _bmat = lambda *arg: np.bmat(arg[-1])
    _diag = lambda *arg: np.diag(arg[-1])
    _eye = lambda *arg: np.eye(arg[-1])

    """Docstring for Class"""
    def __init__(self,T11,T12,T21,T22):
        """Transfer Matrix object
        """
        self.T11,self.T12,self.T21,self.T22 = T11,T12,T21,T22

    def transfer(self,c_in_left,c_out_left):
        """Acts the transfer matrix on a mode coef vector
        """
        return np.asarray(self.T11 @ c_in_left + self.T12 @ c_out_left).flatten(), \
                np.asarray(self.T21 @ c_in_left + self.T22 @ c_out_left).flatten()

    def calculate_SM(self):
        """Calculates the associated scattering matrix
        """
        T22inv = self._inv(self.T22)
        self.S11 = - T22inv @ self.T21
        self.S12 = T22inv
        self.S21 = self.T11 - self.T12 @ T22inv @ self.T21
        self.S22 = self.T12 @ T22inv

        return self

    def scatter(self, c_in_left, c_in_right):
        """Acts the transfer matrix on inbound modes
        """
        try:
            return np.asarray(self.S11 @ c_in_left + self.S12 @ c_in_right).flatten(), \
                    np.asarray(self.S21 @ c_in_left + self.S22 @ c_in_right).flatten()
        except AttributeError as E:
            raise AttributeError("Scattering matrix uncalculated, must call .calculate_SM() first") from E

    def __matmul__(self,other):
        """Combination operator for transfer matricies
        """
        if other.__class__ == TransferMatrix:
            cls = TransferMatrix
        else:
            cls = self.__class__
        return cls(
                    T11 = other.T11 @ self.T11 + other.T12 @ self.T21,
                    T12 = other.T11 @ self.T12 + other.T12 @ self.T22,
                    T21 = other.T21 @ self.T11 + other.T22 @ self.T21,
                    T22 = other.T21 @ self.T12 + other.T22 @ self.T22
            )

    def copy(self):
        """Returns a hard copy of the transfer matrix
        """
        return self.__class__(self.T11.copy(),self.T12.copy(),self.T21.copy(),self.T22.copy())

    @classmethod
    def for_symmetric_layer(cls,mode,mode_gap,L):
        alpha = mode.W_inv @ mode_gap.W
        beta = mode.V_inv @ mode_gap.V
        A = (alpha + beta)/2
        B = (alpha - beta)/2
        alpha = mode_gap.W_inv @ mode.W
        beta = mode_gap.V_inv @ mode.V
        Ap = (alpha + beta)/2
        Bp = (alpha - beta)/2
        Lk0 = 2*np.pi*L/np.append(mode.k.wls,mode.k.wls)
        X = cls._diag(np.exp(mode.lam*Lk0))
        Xs = cls._diag(np.exp(-mode.lam*Lk0))

        XA = X @ A
        XB = X @ B
        XsB = Xs @ B
        XsA = Xs @ A

        return cls(
            Ap @ XA + Bp @ XsB,
            Ap @ XB + Bp @ XsA,
            Bp @ XA + Ap @ XsB,
            Bp @ XB + Ap @ XsA
        )
        # A = mode.W_inv @ mode_gap.W + mode.V_inv @ mode_gap.V
        # A_inv = cls._inv(A)
        # B = mode.W_inv @ mode_gap.W - mode.V_inv @ mode_gap.V
        # Lk0 = 2*np.pi*L/np.append(mode.k.wls,mode.k.wls)
        # X = cls._diag(np.exp(mode.lam*Lk0))
        #
        # G = cls._inv(A - X @ B @ A_inv @ X @ B)
        #
        # S11 = G @ (X @ B @ A_inv @ X @ A - B)
        # S12 = G @ X @ (A - B @ A_inv @ B)
        # S21 = S12
        # S22 = S11
        #
        # S12inv = cls._inv(S12)
        # return cls(
        #         S21 - S22 @ S12 @ S11,
        #         S22 @ S12inv,
        #         - S12inv @ S11,
        #         S12inv
        #     )

    @classmethod
    def for_interface(cls,modeL,modeR):
        A = 0.5 * (modeR.W_inv @ modeL.W + modeR.V_inv @ modeL.V)
        B = 0.5 * (modeR.W_inv @ modeL.W - modeR.V_inv @ modeL.V)

        return cls(A,B,B,A)

    @classmethod
    def for_propegation(cls,lam,zp):
        return cls(
            cls._diag(np.exp(lam*zp)),
            cls._sqzeros(len(lam)),
            cls._sqzeros(len(lam)),
            cls._diag(np.exp(-lam*zp)))

    @classmethod
    def for_mode_to_field(cls,mode):
        return cls(mode.W,mode.W,mode.V,-mode.V)
    @classmethod
    def for_field_to_mode(cls,mode):
        return cls(
            mode.W,
            cls._sqzeros(mode.W.shape[0]),
            cls._sqzeros(mode.W.shape[0]),
            mode.W)

    @classmethod
    def for_null(cls,num_modes=2):
        return cls(
            cls._eye(num_modes),
            cls._sqzeros(num_modes),
            cls._sqzeros(num_modes),
            cls._eye(num_modes)
        )

class TransferMatrixSparse(TransferMatrix):
    _inv = lambda *arg: sla.inv(arg[-1])
    _sqzeros = lambda *arg: 0*sp.eye(arg[-1],format="csc")
    _bmat = lambda *arg: sp.bmat(arg[-1],format="csc")
    _diag = lambda *arg: sp.diags(arg[-1],format="csc")
    _eye = lambda *arg: sp.eye(arg[-1],format="csc")
