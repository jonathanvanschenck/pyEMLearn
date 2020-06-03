"""Docstring for Module"""

from cmath import sqrt as csqrt
from numpy import interp,sqrt

from pyEMLearn2.utils import BiaxialMatrix,UniaxialMatrix

class Material:
    """Docstring for Class"""
    def __init__(self,name=None,isotropic=True):
        """Base Material class, intended to be inhereted from
        
        Parameters
        ----------
        name : str or None, optional
            A name to describe the of material. Default is None.
            
        isotropic : bool
            Boolean indicating if the material is isotropic.
            
        Returns
        -------
        None
        """
        self.name = name
        self.isotropic = isotropic
        
    @property
    def name(self):
        if self._name is None:
            return self.__class__.__name__
        return self._name
    @name.setter
    def name(self,name):
        assert type(name) == str or name is None
        self._name = name
    
    def __repr__(self):
#        if self.name is None:
#            name = ""
#        else:
#            name = ": " + self.name
        return "<Material {}>".format(self.name)
    
    def n(self,wl):
        """Get the complex index of refraction at a specified wavelength

        Parameters
        ----------
        wl : float
            A wavelength value for which to return the index of refraction

        Raises
        ------
        NotImplementedError
            In the base class, this must be overridden

        Returns
        -------
        None.

        """
        raise NotImplementedError
        
    def er(self,wl):
        """Get the complex relative permittivity at a specified wavelength

        Parameters
        ----------
        wl : float
            A wavelength value for which to return the index of refraction

        Returns
        -------
        er : complex
            relative permeability

        """
        n,k = self.n(wl)
        return (n + 1j*k)**2
    
    def mr(self,wl):
        """Get the (complex) relative permeability at a specified wavelength

        Parameters
        ----------
        wl : float
            A wavelength value for which to return the relative permeability

        Returns
        -------
        mr : complex
            relative permeability

        """
        return 1.0+0j
    
    def __add__(self,other):
        return SumMaterial(self,other)
    
class SumMaterial(Material):
    def __init__(self,matA,matB,name=None):
        self.matA = matA
        self.matB = matB
        Material.__init__(self,name=name)
        
    def __repr__(self):
        return "<Material Sum: {0} + {1}>".format(self.matA.name,self.matB.name)
        
    def n(self,wl):
        n = csqrt(self.er(wl))
        return n.real,n.imag
    
    def er(self,wl):
        return self.matA.er(wl) + self.matB.er(wl)

class ConstantIndex(Material):
    """Docstring for Class"""
    def __init__(self,n0,k0=0,name=None):
        """Material with a index constant over wavelength

        Parameters
        ----------
        n0 : float
            real part of the index of refraction
        k0 : float, optional
            imaginary part of the index of refraction. The default is 0.
        name : str or None, optional
            A name to describe the of material. Default is None.
            
        Returns
        -------
        None
        """
        Material.__init__(self,name)
        self.n0 = n0
        self.k0 = k0
        
    def __repr__(self):
        return "<Material {0}: Constant Index: n0={1}, k0={2}>".format(
                self.name,
                self.n0,
                self.k0
            )
        
    def n(self,wl):
        """Get the (complex) index of refraction at a specified wavelength

        Parameters
        ----------
        wl : float
            A wavelength value for which to return the index of refraction

        Returns
        -------
        n_re
            Real part of the index of refraction
        n_im
            Imaginary part of the index of refraction
        """
        return self.n0,self.k0
    
    
class SellmeierIndex(Material):
    """Docstring for Class"""
    def __init__(self,Bcoef,Ccoef,name=None):
        """Material with an index modelled with the Sellmeier equation:
            n(wl) = 1 + SUM_i { Bi / (Ci - wl^-2) }

        Parameters
        ----------
        Bcoef : interable of floats
            List of B coefficents used in Sellmeier equation. Must be the 
            same length as Ccoef
        Ccoef : interable of floats
            List of C coefficents used in Sellmeier equation. Must be the 
            same length as Bcoef
        name : str or None, optional
            A name to describe the of material. Default is None.
            
        Returns
        -------
        None
        """
        assert len(Bcoef)==len(Ccoef), "Coefficent Lists must have same length"
        Material.__init__(self,name=name)
        self.B = Bcoef
        self.C = Ccoef
        
    def __repr__(self):
        return "<Material {0}: Sellmeier: B= {1}, C= {2}>".format(
                  self.name,
                  ", ".join([str(_B) for _B in self.B]),
                  ", ".join([str(_C) for _C in self.C])
              )
        
    def n(self,wl):
        """Get the (complex) index of refraction at a specified wavelength

        Parameters
        ----------
        wl : float
            A wavelength value for which to return the index of refraction

        Returns
        -------
        n_re
            Real part of the index of refraction
        n_im
            Imaginary part of the index of refraction
        """
        n_re = 1
        for _B,_C in zip(self.B,self.C):
            n_re += _B/(_C-wl**(-2.0))
        return n_re,0.0
    
class SellmeierIndexSquare(Material):
    """Docstring for Class"""
    def __init__(self,Bcoef,Ccoef,name=None):
        """Material with an index modelled with the "square" Sellmeier equation:
            n(wl)^2 = 1 + SUM_i { Bi wl^2 / (wl^2 - Ci) }

        Parameters
        ----------
        Bcoef : interable of floats
            List of B coefficents used in Sellmeier equation. Must be the 
            same length as Ccoef
        Ccoef : interable of floats
            List of C coefficents used in Sellmeier equation. Must be the 
            same length as Bcoef
        name : str or None, optional
            A name to describe the of material. Default is None.
            
        Returns
        -------
        None
        """
        Material.__init__(self,name)
        self.B = Bcoef
        self.C = Ccoef
        assert len(self.B)==len(self.C), "Coefficent Lists must have same length"
        
    def __repr__(self):
        return "<Material {0}: 'square' Sellmeier: B= {1}, C= {2}>".format(
                  self.name,
                  ", ".join([str(_B) for _B in self.B]),
                  ", ".join([str(_C) for _C in self.C])
              )
        
    def n(self,wl):
        """Get the (complex) index of refraction at a specified wavelength

        Parameters
        ----------
        wl : float
            A wavelength value for which to return the index of refraction

        Returns
        -------
        n_re
            Real part of the index of refraction
        n_im
            Imaginary part of the index of refraction
        """
        n_resq = 1
        for _B,_C in zip(self.B,self.C):
            n_resq += _B*wl**2/(wl**2-_C)
        return sqrt(n_resq),0.0
    
class CauchyIndex(Material):
    """Docstring for Class"""
    def __init__(self,coef,name=None):
        """Material with an index modelled with the Cauchy equation:
            n(wl) = c0 + c1 / wl^2 + c2 / wl^4 + ...  

        Parameters
        ----------
        coef : interable of floats
            List of coefficients used in Cauchy equation.
        name : str or None, optional
            A name to describe the of material. Default is None.
            
        Returns
        -------
        None
        """
        Material.__init__(self,name)
        self.coef = coef
        
    def __repr__(self):
        return "<Material {0}: Cauchy: {1}>".format(
                  self.name,
                  ", ".join(["c_{0} = {1}".format(i,_c) for i,_c in enumerate(self.coef)])
              )
        
    def n(self,wl):
        """Get the (complex) index of refraction at a specified wavelength

        Parameters
        ----------
        wl : float
            A wavelength value for which to return the index of refraction

        Returns
        -------
        n_re
            Real part of the index of refraction
        n_im
            Imaginary part of the index of refraction
        """
        n_re = 0
        for i,c in enumerate(self.coef):
            n_re += c/wl**(2.0*i)
        return n_re,0.0
    
class DiscreteIndex(Material):
    """Docstring for Class"""
    def __init__(self,wls_n,wls_k=None,name=None):
        """Material with an index modelled interpolating from a provided table

        Parameters
        ----------
        wls_n : numpy array
            A numpy array containing the wavelengths and real parts of the 
            index of refraction. Array must be size (m x 2), with the first
            column containing wavelength values, and the second column 
            containing index values. Need not be the same size as wls_k, but
            must have the same wavelength units as in wls_k
        wls_k : numpy array, optional
            A numpy array containing the wavelengths and imaginary parts of the 
            index of refraction. Array must be size (m x 2), with the first
            column containing wavelength values, and the second column 
            containing index values. Need not be the same size as wls_n, but
            must have the same wavelength units as in wls_n. If ommited, the
            default imaginary part is set to 0.
        name : str or None, optional
            A name to describe the of material. Default is None.
            
        Returns
        -------
        None
        """
        Material.__init__(self,name=name)
        self.n_re = lambda wl: interp(wl, wls_n[:,0],wls_n[:,1])
        if wls_k is None:
            self.n_im = lambda wl: 0.0
        else:
            self.n_im = lambda wl: interp(wl, wls_k[:,0],wls_k[:,1])
    
    def __repr__(self):
        return "<Material {0}: Discrete: n_re={1}, n_im={2}>".format(
                  self.name,
                  self.n_re,
                  self.n_im
              )
    
    def n(self,wl):
        """Get the (complex) index of refraction at a specified wavelength

        Parameters
        ----------
        wl : float
            A wavelength value for which to return the index of refraction

        Returns
        -------
        n_re
            Real part of the index of refraction
        n_im
            Imaginary part of the index of refraction
        """
        return self.n_re(wl),self.n_im(wl)

class LorentzIndex(Material):
    """Docstring for Class"""
    def __init__(self,A=[],G=[],E=[],eps_inf=1.0,name=None):
        """Material with an index modelled with the Drude-Lorentz equation:
            e_r(E) = eps_inf + SUM_i { Ai / (Ei^2 - Gi*E*j - E^2) }
            n(E) = sqrt(e_r)

        Parameters
        ----------
        A : interable of floats, optional
            List of A values (controls oscillator strength) used in DL 
            equation. Must be the same length as G and E. The default is [].
        G : interable of floats, optional
            List of G values (controls broadening) in eV used in DL 
            equation. Must be the same length as A and E. The default is [].
        E : interable of floats, optional
            List of E values (controls center energy of transition) in eV used
            in DL equation. Must be the same length as G and A. 
            The default is [].
        eps_inf : float, optional
            Value of the dielectric constant at infinity. The default is 1.0.
        name : str or None, optional
            A name to describe the of material. Default is None.
            
        Returns
        -------
        None
        """
        Material.__init__(self,name=name)
        assert len(A)==len(G), "Coefficent Lists must have same length"
        assert len(A)==len(E), "Coefficent Lists must have same length"
        self.A = A
        self.G = G
        self.E = E
        self.einf = eps_inf
    
    def __repr__(self):
        return "<Material {0}: Lorentz: eps_inf={1}; {2}>".format(
                self.name,
                self.einf,
                "; ".join(["A = {0}, G = {1} eV, E = {2} eV".format(_A,_G,_E)\
                             for _A,_G,_E in zip(self.A,self.G,self.E)])
              )
        
    def n(self,wl):
        """Get the (complex) index of refraction at a specified wavelength (in um)

        Parameters
        ----------
        wl : float
            A wavelength value (in um) for which to return the index of refraction

        Returns
        -------
        n_re
            Real part of the index of refraction
        n_im
            Imaginary part of the index of refraction
        """
        e = 1.238/wl
        er = 1*self.einf
        for _A, _G, _E in zip(self.A,self.G,self.E):
            er += _A/(_E**2 - 1j*_G*e - e**2)
        n = csqrt(er)
        return n.real,n.imag
    
class LorentzIndex2(LorentzIndex):
    """Docstring for Class"""
    def __init__(self,A=[],G=[],E=[],eps_inf=1.0,name=None):
        """Material with an index modelled with the Drude-Lorentz equation:
            e_r(E) = eps_inf + SUM_i { Ai*Ei*Gi / (Ei^2 - Gi*E*j - E^2) }
            n(E) = sqrt(e_r)

        Parameters
        ----------
        A : interable of floats, optional
            List of A values (controls oscillator strength) used in DL 
            equation. Must be the same length as G and E. The default is [].
        G : interable of floats, optional
            List of G values (controls broadening) in eV used in DL 
            equation. Must be the same length as A and E. The default is [].
        E : interable of floats, optional
            List of E values (controls center energy of transition) in eV used
            in DL equation. Must be the same length as G and A. 
            The default is [].
        eps_inf : float, optional
            Value of the dielectric constant at infinity. The default is 1.0.
        name : str or None, optional
            A name to describe the of material. Default is None.
            
        Returns
        -------
        None
        """
        LorentzIndex.__init__(self,A=A,G=G,E=E,eps_inf=eps_inf,name=name)
        
    def n(self,wl):
        e = 1.238/wl
        er = 1*self.einf
        for _A, _G, _E in zip(self.A,self.G,self.E):
            er += _A*_E*_G/(_E**2 - 1j*_G*e - e**2)
        n = csqrt(er)
        return n.real,n.imag

class BiaxialMaterial(Material):
    def __init__(
        self,
        mat_x,mat_y,mat_z,
        ex,ey,ez,
        name=None,
    ):
        """docstring"""
        Material.__init__(self,name,isotropic=False)
        self.materials = [mat_x,mat_y,mat_z]
        self.ERMatrix = BiaxialMatrix(ex,ey,ez)
        
    def er(self,wl):
        return self.ERMatrix(*[mat.er(wl) for mat in self.materials])
    def mr(self,wl):
        return np.diag([1.,1.,1.])+0j