"""
Full Cosmology class
Includes old PowSpec class functions, grouped together for easier use
"""

"""
(MILDLY) USEFUL COSMOLOGY CLASS.
NOTE ALL UNITS ARE WRT h^-1

Much of the cosmography is from (as ever) Hogg arXiv:astro-ph/9905116

(C) IAN HARRISON
2012-
IAN.HARRISON@ASTRO.CF.AC.UK

"""



import numpy as np
from numpy import sqrt, log, exp, fabs, pi, cos, sin

from scipy import integrate
from scipy import interpolate
from matplotlib import pyplot as plt

from constants import *
import hmf
import pickle_script
import powspec as ps


class Cosmology:

  # default gives a dictionary of default cosmology initialisation values
  def __init__(self, ):
    """
    Constructor.
    Default initialisation with WMAP7+BAO+H0 parameters,
    Tinker mass function and fitted P(k).
    """

    default={"dc":1.686e0, "h0":0.702e0, "om":0.274e0, "ode":0.725e0,
         "w0":-1.e0, "ob":0.0458, "o_r":8.6e-5, "o_k":0.e0, "f_nl":0.e0,
         "tau_r":0.087, "z_r":10.4, "ns":0.96,
         "hmf_type":'tinker'}

    # create / load pickle of given cosmology
    pickle_script.init_pickle(self,default)

    print "Cosmology \n ......... \n"

    self.eta_r=self.eta(self.z_rec)

    #self.mf = hmf.Hmf(mf_type=hmf_type)

  def display(self):
    """
    Displays what you're working with.
    """
    print("YOU ARE USING THE FOLLOWING COSMOLOGY:")
    print("{h_0, O_m, O_de, w_0, O_b, O_r, O_k} = ")
    print("{{{}, {}, {}, {}, {}, {}, {}}}".format(self.h_0, self.O_m0, self.O_de0, self.w_0, self.O_b0, self.O_r0, self.O_k0))
    print("WITH MASS FUNCTION:")
    self.mf.display()
    print("WITH POWER SPECTRUM:")
    self.pk.display()

  def set_hmf(self, set_mf):
    """
    Set method for the HMF within the Cosmology
    """
    self.mf = set_mf

  def set_powspec(self, set_pk):
    """
    Set method for power spectrum within the Cosmology
    """
    self.pk = set_pk

  def rho_c(self):
    """
    Critical density
    """
    return rhofactor * ((3.e0*100.e0*100.e0) / (8.e0*pi*G*Hfactor*Hfactor))

  def rho_m(self, z):
    """
    Average density of Universe at redshift z
    """
    return (rhofactor*self.O_m(z) *
            ((3.e0*100.e0*100.e0) / (8.e0*pi*G*Hfactor*Hfactor)))

  def O_m(self, z):
    """
    Omega matter.
    """
    return (self.O_m0*pow(1.e0+z, 3.e0) /
           (self.O_m0*pow(1.e0+z, 3.e0) + 1.e0 - self.O_m0));

  def h(self, z):
    """
    Hubble function ***little*** h(z).
    """
    a = 1.e0/(1.e0+z)
    h_square = (self.O_r0/pow(a, 4.e0) +
                self.O_m0/pow(a, 3.e0) +
                self.O_k0/pow(a, 2.e0) +
                self.O_de0/pow(a, 3.e0*(1.e0+self.w_0)))

    return sqrt(h_square)

  def _integ_one_on_h_scalar(self, z_min, z_max):
    """Scalar calculation of 1/h between z_min and z_max
    """
    integfunc = lambda x : 1.e0/self.h(x)
    return integrate.quad(integfunc, z_min, z_max)[0]

  _integ_one_on_h_vector = np.vectorize(_integ_one_on_h_scalar)

  def H(self, z):
    """
    Hubble function ***upper*** H(z).
    """
    return self.H_0*self.h(z)

  def dvdz(self, z):
    """
    Comoving volume element at redshift z.
    """
    return 4.e0*2998e0*(1+z)*(1+z)*pi*pow(self.D_a(z),2.e0) / self.h(z)

  def V_between(self, z_min, z_max):
    """
    Volume between two redshifts.
    """
    points = 200
    int_arr = np.zeros(points)
    z_arr = np.linspace(z_min, z_max, points)
    for i in np.arange(points):
      int_arr[i] = self.dvdz(z_arr[i])

    return integrate.trapz(int_arr, z_arr)

  def D_a(self, z):
    """
    Angular diameter distance to redshift z.
    """
    if np.isscalar(z):
      integ = self._integ_one_on_h_scalar(0.e0, z)
    else:
      integ = self._integ_one_on_h_vector(self, 0.e0, z)

    return 2998.e0*integ/(1.e0 + z)

  def D_l(self, z):
    """
    Luminosity distance to redshift z.
    """
    return D_a(z)*(1.e0 + z)**2.e0

  def D_c(self, z):
    """
    Comoving radial distance to redshift z.
    """
    return D_a(z)*(1.e0 + z)/(100.e0*self.h0)

  def dist_mod(self, z):
    """
    Distance modulus.
    5 * log(D_l(z)) + 25
    """
    return 5.e0*np.log10(self.D_l(z)) + 25.e0

  def eta(self, z):
    """
    Size of particle horizon at redshift z
    """
    if (z == np.inf):
      retVar = 0.e0
    else:
     if np.isscalar(z):
        retVar = self._integ_one_on_h_scalar(z, np.inf)/(100.e0*self.h_0)
     else:
        retVar = self._integ_one_on_h_vector(self, z, np.inf)/(100.e0*self.h_0)
    return retVar

  def _integ_lookback_scalar(self, z):
    """Scalar calculation of lookback time to redshift z
    """
    integfunc = lambda x : 1.e0/(self.h(x)*(1.e0 + x))
    return integrate.quad(integfunc, 0.e0, z)[0]/(100.e0*self.h_0)

  _integ_lookback_vector = np.vectorize(_integ_lookback_scalar)

  def t_lookback(self, z):
    """Lookback time to redshift z
    """
    if np.isscalar(z):
      t_lb = self._integ_lookback_scalar(z)
    else:
      t_lb = self._integ_lookback_vector(self, z)
    return t_lb

  def _integ_growth_scalar(self, z):
    """Scalar calculation of linear growth function
    """
    integfunc = lambda x : 1.e0/pow(self.h((1.e0/x) -1.e0)*x, 3.e0)

    return ((5.e0/2.e0) * self.O_m0 * self.h(z) *
            integrate.quad(integfunc, 0.e0, 1.e0/(1.e0+z))[0])

  _integ_growth_vector = np.vectorize(_integ_growth_scalar)

  def growth(self, z):
    """Linear growth function.
    """
    if np.isscalar(z):
      self.Dplus = self._integ_growth_scalar(z)
    else:
      self.Dplus = self._integ_growth_vector(self, z)
    return self.Dplus


  def dndlnm(self, lnm, z):
    """
    Comoving number density of dark matter haloes in logarithmic m.
    """
    s = self.sigma_wmap7fit(lnm)
    dlnsdlnm = self.dlnsigmadlnm_wmap7fit(lnm)
    D_z = self.growth(z)
    D_0 = self.growth(0.e0)

    return (self.mf.f_sigma(s*D_z/D_0, z) *
            fabs(log(D_z/D_0) + dlnsdlnm) *
            self.rho_m(z) / (exp(lnm)) *
            self.mf.r_ng(s*D_z/D_0, self.delta_c, self.fnl))

  def computeNinBin(self, z_min, z_max, lnm_min, lnm_max):
    """
    Total number of dark matter haloes expected within a given mass,
    redshift bin.
    """
    return integrate.dblquad(self.dNdlnmdz,
                             z_min, z_max,
                             lambda lnm: lnm_min, lambda lnm: lnm_max)[0]

  def computeLittleNinZBin(self, lnm_min, lnm_max, z):
    """
    Total number of haloes within a given mass bin at fixed redshift.
    """
    return integrate.quad(self.dndlnm, lnm_min, lnm_max, args=(z))[0]

  def dNdlnmdz(self, lnm, z):
    """
    Total number of dark matter haloes at a given redshift.
    Product of dndlnm * dVdz
    """
    return self.dndlnm(lnm, z)*self.dvdz(z)

  def dNdlnm0dz(self, lnm, z):
    """
    Total number of dark matter haloes,
    of equivalent mass at redshift zero, at a given redshift.
    Product of dndlnm * dVdz
    """
    return self.dndlnm(lnm, 0.e0)*self.dvdz(z)


















  def display(self):
    """Display method to show power spectrum currently working with.
    """
    print("Power spectrum {}".format(self.label))


  def vd_initialisation(self,rrange,mrange,z=0.0):
    """ initialise parameters required for
        void_distribution.py script """

    self.growth(z)
    self.sigma_r(rrange)
    self.sig_fit_z = self.sigma_fit(rrange,self.sigmar_z)
    self.sig_fit = self.sigma_fit(rrange,self.sigmar)
    self.dlnsigma_dlnr(rrange,z)
    self.dlnsigma_dlnm(mrange,z)

    return None

  def sigma_fit(self,rrange,sigma_r):
    """ fit to sigma """
    fit = np.polyfit(log(rrange),sigma_r,6)
    sig_fit = np.poly1d(fit)

    return sig_fit

  def sigma_r(self,r):
    """ returns root of the matter variance, smoothed
        with a top hat window function at a radius r
        NOTE self.sigmar is NOT redshift dependent
        self.sigmar_z IS redshift dependent """

    if np.isscalar(r):
      s_r_sq, s_r_sq_error = self.sigma_r_sq(r)
      self.sigmar_z = sqrt(s_r_sq) * self.Dplus
      self.sigmar = sqrt(s_r_sq)
    else:
      s_r_sq, s_r_sq_error = self.sigma_r_sq_vec(self,r)
      self.sigmar_z = sqrt(s_r_sq) * self.Dplus
      self.sigmar = sqrt(s_r_sq)

    self.sigmar.tolist()
    return self.sigmar

  def sigma_r_sq(self,r):
    """ integrate the function in sigma_integral
        between the limits of k : 0 to inf. """

    s_r_sq, s_r_sq_error = integrate.quad(self.sigma_integral,0.,np.inf,args=(r))#,limit=10000)

    return s_r_sq, s_r_sq_error

  sigma_r_sq_vec = np.vectorize(sigma_r_sq)    #vectorize sigma-squared function

  def sigma_integral(self,k,r,z=0.0):
    """ returns the integral required to calculate sigma
        squared (Coles & Lucchin pg.266, A.Zentner 06 eq.14)"""

    return (k**2 / (2 * pi**2)) * fabs(self.tophat_w(k,r)) * ps.power_spectrum_P(self,k,z)

  def tophat_w(self, k, r):
    """ Fourier transform of the real space tophat
        window function (eq.9 from A.Zentner 06) """

    return (3.*(sin(k*r) - k*r*cos(k*r)))/((k*r)**3.)


  def dlnsigma_dlnr(self,rrange):
    """ slope of root matter variance wrt log radius:
        d(log(sigma)) / d(log(r))

        Polynomial fit to supplied cosmology.
        Returns poly1d object """

    self.sigma_r(rrange)

    fit = np.polyfit(log(rrange),log(self.sigmar_z),6)
    self.dlnsigmadlnr = np.polyder(np.poly1d(fit))

    return self.dlnsigmadlnr


  def dlnsigma_dlnm(self,mrange,z):
    """ slope of root matter variance wrt log mass:
        d(log(sigma)) / d(log(M))

        Polynomial fit to supplied cosmology.
        Returns poly1d object """

    rho = self.rho_m(z)

    # convert mass -> radius to get sigma
    rrange = ((3*mrange)/(4*pi*rho))**(0.333333333333333333333)
    self.sigma_r(rrange)

    fit = np.polyfit(log(mrange),log(self.sigmar_z),6)
    self.dlnsigmadlnm = np.polyder(np.poly1d(fit))

    return self.dlnsigmadlnm


  def sigma_wmap7fit(self, lnm):
    """
    Root of matter variance smoothed with top hat window function on a scale
    specified by log(m)

    Polynomial fit to calculation from a CAMB power spectrum
    with WMAP7 parameters
    """
    return np.exp(18.0058 - 1.47523*lnm + 0.0385449*lnm*lnm - 0.0000112539*pow(lnm,4) + (1.3719e-9)*pow(lnm,6))


  def dlnsigmadlnm_wmap7fit(self, lnm):
    """
    Slope of root matter variance wrt log mass:
    d(log(sigma)) / d(log(m))

    Polynomial fit to calculation from a CAMB power spectrum
    with WMAP7 parameters
    """
    return -1.47523 + 0.0770898*lnm - 0.0000450156*pow(lnm,3) + (8.23139e-9)*pow(lnm,5)



