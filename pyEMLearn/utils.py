"""pyEMLearn utilities module

This module contains various objects useful for other modules, like quaterion-
based rotations, Transfer/Scattering Matricies classes, etc.
"""

from numpy import eye,diag,exp,zeros,array,append,\
                  sqrt,dot,cross,sin,cos,pi,arccos,\
                  arcsin, transpose
from numpy import all as _all
from numpy.linalg import inv, norm, det

class Quaternion:
    def __init__(self,real,imag):
        self.real = real
        assert len(imag) == 3
        self.imag = array(imag)

    def __repr__(self):
        if _all(self.imag == 0.0):
            s = "{:+}".format(self.real)
        elif float(self.real) == 0.0:
            s = "{0}i{1}j{2}k".format(*["{:+}".format(e) for e in self.imag])
        else:
            s = "{0}{1}i{2}j{3}k".format(*["{:+}".format(e) for e in self.as_array()])
        if s[0] == "+":
            s = s[1:]
        return s


    def as_array(self):
        return append([self.real],self.imag)

    def __eq__(self,other):
        if type(other) == Quaternion:
            return (self.real==other.real and self.imag==other.image)
        else:
            return False

    def __ne__(self,other):
        return not self.__eq__(other)

    def __add__(self,other):
        if type(other) == Quaternion:
            return Quaternion(self.real+other.real,
                             self.imag+other.image)
        else:
            return Quaternion(self.real+other,
                             self.imag)
    def __radd__(self,other):
        # Only here if other isn't a quaterion
        return Quaternion(other+self.real,
                         self.imag)

    def __sub__(self,other):
        if type(other) == Quaternion:
            return Quaternion(self.real-other.real,
                             self.imag-other.image)
        else:
            return Quaternion(self.real-other,
                             self.imag)
    def __rsub__(self,other):
        # Only here if other isn't a quaterion
        return Quaternion(other-self.real,
                         self.imag)

    def __neg__(self):
        return Quaternion(-self.real,-self.imag)
    def __pos__(self):
        return Quaternion(1*self.real,self.imag.copy())
    def conj(self):
        return Quaternion(self.real,-self.imag)
    def __invert__(self):
        return self.conj()
    def __abs__(self):
        return sqrt(self.real*self.real + dot(self.imag,self.imag))

    def __mul__(self,other):
        if type(other) == Quaternion:
            return Quaternion(
                self.real*other.real - dot(self.imag,other.imag),
                self.real*other.imag + other.real*self.imag + cross(self.imag,other.imag)
            )
        else:
            return Quaternion(other*self.real,other*self.imag)

    def __rmul__(self,other):
        # Only here if other isn't a quaterion
        return Quaternion(other*self.real,other*self.imag)

    def __truediv__(self,other):
        if type(other) == Quaternion:
            raise TypeError
        else:
            return Quaternion(self.real/other,self.imag/other)

    def __getitem__(self, key):
        try:
            ik = int(key)
        except:
            raise KeyError
        if ik < -4:
            raise IndexError
        elif ik == 0 or ik == -4:
            return self.real
        elif ik < 0:
            return self.imag[ik]
        else:
            return self.imag[ik-1]

    def __setitem__(self, key, value):
        try:
            fval = float(value)
        except:
            raise ValueError
        try:
            ik = int(key)
        except:
            raise KeyError
        if ik < -4:
            raise IndexError
        elif ik == 0 or ik == -4:
            self.real = fval
        elif ik < 0:
            self.imag[ik] = fval
        else:
            self.imag[ik-1] = fval

    @classmethod
    def from_vector(cls,vector):
        assert len(vector) == 4
        return cls(vector[0],vector[1:])


class Rotation3D:
    def __init__(self,angle=0,axis=[0,0,1]):
        l = norm(axis)
        assert l > 0
        naxis = array(axis)/l
        s = sin(angle*pi/360)
        c = cos(angle*pi/360)
        self.q = Quaternion(c,s*naxis)
        self.q_inv = ~self.q
        self.v = Quaternion(0,[0,0,0])

    def get_angle(self):
        return 360 * arccos(self.q.real) / pi
    def get_axis(self):
        s = sin(arccos(self.q.real))
        return self.q.imag/s

    @classmethod
    def from_q(cls,q):
        assert type(q) == Quaternion
        assert abs(abs(q)-1) < 1e-6
        newR3D = cls()
        newR3D.q = 1*q
        newR3D.q_inv = ~newR3D.q
        return newR3D

    @classmethod
    def from_euler(cls,alpha,beta,gamma):
        """x,y,z are external coords, XYZ are rotated coords, N is intersection of Z* with xy
        alpha gives angle between N and x (precession)
        beta gives angle between Z and z (nutation)
        gamma gives angle between X and N (intrinsic rotation)
        """
        return cls(alpha,[0,0,1]) @ cls(beta,[1,0,0]) @ cls(gamma,[0,0,1])

    def __call__(self,vec):
        assert len(vec) == 3
        self.v.imag = array(vec)
        return (self.q * self.v * self.q_inv).imag

    def __matmul__(self,other):
        if type(other) == Rotation3D:
            return Rotation3D.from_q(self.q * other.q)
        else:
            raise TypeError


# Not currently used
class FreezeableParameters:
    def __init__(self,p,pname):
        assert len(p) == len(pname)
        self.p = array(p)
        self.pname = array(pname)
        self.freeze_all()
        self.pnum = len(self.p)

    def __repr__(self):
        return "<FreezeableParameters: Active = {}>".format(", ".join(self.pname[self.which]))

    def freeze_all(self):
        self.which = array(len(self.pnum)*[False])

    def modify_parameter(self,pname,active=True):
        if type(pname) == str:
            try:
                i = list(self.pname).index(pname)
            except ValueError:
                raise ValueError("Invalid parameter name")
            self.which[i] = active
            return None
        try:
            i = int(pname)
        except ValueError:
            raise TypeError("Invalid parameter type")
        self.which[i] = active
        return None

    def get_active_parameters(self):
        return self.p[self.which]

    def get_active_parameter_names(self):
        return self.pname[self.which]

    def update(self,pp):
        assert len(pp) == len(self)
        self.p[self.which] = array(pp)


class BiaxialMatrix:
    def __init__(self,ex,ey,ez):
        """Eigenvectors for biaxial matrix"""
        self.W_inv = array([ex,ey,ez])
        self.W = transpose(self.W_inv)

    def __call__(self,lamX,lamY,lamZ):
        """eigenvalues for biaxial matrix"""
        return self.W @ diag(array([lamX,lamY,lamZ])) @ self.W_inv

    @classmethod
    def from_euler(cls,alpha,beta,gamma):
        r = Rotation3D.from_euler(alpha,beta,gamma)
        return cls(*[r(vec) for vec in [[1,0,0],[0,1,0],[0,0,1]]])

    @classmethod
    def from_rotation(cls,rotation):
        return cls(*[rotation(vec) for vec in [[1,0,0],[0,1,0],[0,0,1]]])

    @classmethod
    def from_angle_axis(cls,angle,axis):
        r = Rotation3D(angle,axis)
        return cls(*[r(vec) for vec in [[1,0,0],[0,1,0],[0,0,1]]])

class UniaxialMatrix(BiaxialMatrix):
    def __init__(self,ex,ey,ez):
        BiaxialMatrix.__init__(self,ex,ey,ez)

    def __call__(self,lamOrd,lamExt):
        """eigenvalues for uniaxial matrix"""
        return self.W @ diag(array([lamExt,lamOrd,lamOrd])) @ self.W_inv

    @classmethod
    def from_extrodinary_axis(cls,extaxis):
        l = norm(extaxis)
        assert l > 0
        rvec = cross([1,0,0],extaxis/l)
        angle = arcsin(norm(rvec))*180/pi
        r = Rotation3D(angle,rvec)
        return cls(*[r(vec) for vec in [[1,0,0],[0,1,0],[0,0,1]]])


class ScatteringMatrix:
    """Docstring for Class"""
    def __init__(self,S11,S12,S21,S22):
        """Scattering Matrix object
        """
        self.S11,self.S12,self.S21,self.S22 = S11,S12,S21,S22

    def scatter(self,c_in_left,c_in_right):
        """Acts the scattering matrix on a mode coef vector
        """
        return self.S11 @ c_in_left + self.S12 @ c_in_right, \
                self.S21 @ c_in_left + self.S22 @ c_in_right


    def __matmul__(self,other):
        """Combination operator scattering matricies (Redheffer Star Product)
        """
        return ScatteringMatrix(
                    self.S11 + self.S12@inv(eye(2)-other.S11@self.S22)@other.S11@self.S21,
                    self.S12@inv(eye(2)-other.S11@self.S22)@other.S12,
                    other.S21@inv(eye(2)-self.S22@other.S11)@self.S21,
                    other.S22 + other.S21@inv(eye(2)-self.S22@other.S11)@self.S22@other.S12
            )

    def copy(self):
        """Returns a hard copy of the scattering matrix
        """
        return ScatteringMatrix(self.S11.copy(),self.S12.copy(),self.S21.copy(),self.S22.copy())

    def det(self):
        """Returns the determinant of the scattering matrix
        """
        return det(vstack([hstack([self.S11,self.S12]),hstack([self.S21,self.S22])]))

    def get_TM(self):
        """Calculates the corresponding transfer matrix
        """
        S12inv = inv(self.S12)
        return TransferMatrix (
                self.S21 - self.S22 * S12inv @ self.S11,
                self.S22 @ S12inv,
                - S12inv @ self.S11,
                S12inv
            )

    @classmethod
    def for_symmetric_layer(cls,W_inv,V_inv,Wg,Vg,lam,k0,L):
        A = W_inv @ Wg + V_inv @ Vg
        A_inv = inv(A)
        B = W_inv @ Wg - V_inv @ Vg
        X = diag(exp(lam*k0*L))
        G = inv(A - X @ B @ A_inv @ X @ B)

        S11 = G @ (X @ B @ A_inv @ X @ A - B)
        S12 = G @ X @ (A - B @ A_inv @ B)

        return cls(S11,S12,S12,S11)

    @classmethod
    def for_interface(cls,WL,VL,WR_inv,VR_inv):
        A = WR_inv @ WL + VR_inv @ VL
        A_inv = inv(A)
        B = WR_inv @ WL - VR_inv @ VL

        return cls(
            B @ A_inv,                 # S11
            0.5 * (A - B @ A_inv @ B), # S12
            2 * A_inv,                 # S21
            - A_inv @ B,               # S22
        )



class TransferMatrix:
    """Docstring for Class"""
    def __init__(self,T11,T12,T21,T22):
        """Transfer Matrix object
        """
        self.T11,self.T12,self.T21,self.T22 = T11,T12,T21,T22

    def transfer(self,c_in_left,c_out_left):
        """Acts the transfer matrix on a mode coef vector
        """

        return self.T11 @ c_in_left + self.T12 @ c_out_left, \
                self.T21 @ c_in_left + self.T22 @ c_out_left


    def __matmul__(self,other):
        """Combination operator for transfer matricies
        """
        return TransferMatrix(
                    T11 = other.T11 @ self.T11 + other.T12 @ self.T21,
                    T12 = other.T11 @ self.T12 + other.T12 @ self.T22,
                    T21 = other.T21 @ self.T11 + other.T22 @ self.T21,
                    T22 = other.T21 @ self.T12 + other.T22 @ self.T22
            )

    def copy(self):
        """Returns a hard copy of the transfer matrix
        """
        return TransferMatrix(self.T11.copy(),self.T12.copy(),self.T21.copy(),self.T22.copy())

    def det(self):
        """Returns the determinant of the transfer matrix
        """
        return det(vstack([hstack([self.T11,self.T12]),hstack([self.T21,self.T22])]))

    def get_SM(self):
        """Calculates the corresponding scattering matrix
        """
        T22inv = inv(self.T22)
        return ScatteringMatrix (
                self.T11 - self.T12 @ T22inv @ self.T21,
                self.T12 @ T22inv,
                - T22inv @ self.T21,
                T22inv
            )
    @classmethod
    def for_symmetric_layer(cls,W_inv,V_inv,Wg,Vg,lam,k0,L):
        raise NotImplementedError


    @classmethod
    def for_interface(cls,WL,VL,WR_inv,VR_inv):
        A = 0.5 * (WR_inv @ WL + VR_inv @ VL)
        B = 0.5 * (WR_inv @ WL - VR_inv @ VL)

        return cls(A,B,B,A)

    @classmethod
    def for_propegation(cls,lam,k0,z):
        return cls(diag(exp(lam*k0*z)),zeros((2,2)),zeros((2,2)),diag(exp(-lam*k0*z)))

    @classmethod
    def for_conversion(cls,W,V):
        return cls(W,W,V,-V)

# class PlaneWaveSource:
#     def __init__(self,kx,ky,kz):
#         self.kx, self.ky, self.kz = kx, ky, kz
#
#
#     def pol_vec(self,AOI,azimuth=0):
#         pass
#
#     @classmethod
#     def from_AOI(cls):
#         pass
#
#     @classmethod
#     def from_AOI(cls):
#         pass
