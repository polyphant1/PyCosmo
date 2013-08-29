"""Void Distribution Script

Python script for reproducing the
number density distribution of voids
with respect to mass and radius.

Also calculates the Extreme Value
Statistics of maximum and minimum
radius voids.
"""

import cPickle as pickle

import math
import numpy as np
from numpy import sqrt, log, exp, fabs, pi, sin, cos
from matplotlib import pyplot as plt
from scipy import integrate
from scipy.interpolate import interp1d

from main import Cosmology
import evs_cosmology as evs

class Void:

  def __init__(self,cosmology):

    self.void_barrier = -2.7
    self.collapse_barrier = 1.686


  def void_and_cloud(self):
    """The Void-and-Cloud parameter

    Parameters
    ----------

    void_barrier : float
    the critical underdensity; defines a void

    collapse_barrier : float
    the critical overdensity; boundary at which virialized objects are formed

    Notes
    -----

    Calculates D, the void-and-cloud parameter, which parametrises
    the impact of halo evolution on the evolving population of voids
    """
    return fabs(self.void_barrier) / \
           (self.collapse_barrier + fabs(self.void_barrier))


  def multiplicity_function_svdw(self,nu,D):
    """Multiplicity function defined in Sheth &
       van de Weyagert (2004 MNRAS 350 517)

    Parameters
    ----------

    nu : array
    scaled void size/mass

    D : float
    void-and-cloud parameter

    void_barrier : float
    the critical underdensity; defines a void

    collapse_barrier : float
    the critical overdensity; boundary at which virialized objects are formed

    Notes
    -----

    calculates equation [4] in Sheth & van de Weygaert
    approximating the infinite series in equation [1]
    """

    return (1/nu) * (nu/(2*pi))**0.5 * exp(-0.5*nu) \
    * exp((-1*fabs(self.void_barrier)*(D**2) \
    / (fabs(self.collapse_barrier)*nu*4)) - (2*(D**4)/(nu**2)))


  multiplicity_function_svdw_vec = np.vectorize(multiplicity_function_svdw)

  def multiplicity_function_jlh(self,sigma,D):
    """Multiplicity function approximation defined
       in Jennings, Li & Hu (2013 MNRAS 000 1)

    Parameters
    ----------

    sigma : array
    standard deviation of the power spectrum

    D : float
    void-and-cloud parameter

    void_barrier : float
    the critical underdensity; defines a void

    collapse_barrier : float
    the critical overdensity; boundary at which virialized objects are formed

    Notes
    -----

    calculates equation [8] in Jennings, Li & Hu
    approximating the infinite series in equation [6]
    """
    x = ((D*sigma)/fabs(self.void_barrier))

    if x <= 0.276:
      return (2./pi)**0.5 * (fabs(self.void_barrier)/sigma) \
             * exp((-1*self.void_barrier**2)/(2.*sigma**2))

    else:
      summation = 0.0
      for j in range(1,5):
        summation = summation + (exp(-0.5*(j*pi*x)**2) *j*pi*x**2.*sin(j*pi*D))

      return 2. * summation

  multiplicity_function_jlh_vec = np.vectorize(multiplicity_function_jlh)

  def multiplicity_function_jlh_exact(self,sigma,D):
    """Exact multiplicity function defined
       in Jennings, Li & Hu (2013 MNRAS 000 1)

    Parameters
    ----------

    sigma : array
    standard deviation of the power spectrum

    D : float
    void-and-cloud parameter

    void_barrier : float
    the critical underdensity; defines a void

    collapse_barrier : float
    the critical overdensity; boundary at which virialized objects are formed

    Notes
    -----

    calculates equation [6] in Jennings, Li & Hu; the exact
    """
    x = ((D*sigma)/fabs(self.void_barrier))

    summation = 0.0
    for j in range(1,500):
      summation = summation + (exp(-0.5*(j*pi*x)**2) *j*pi*x**2.*sin(j*pi*D))

    return 2. * summation


  def multiplicity_function_pls(self,sigma):
    """ multiplicity function for a single crossing barrier from
    Paranjape, Lam & Sheth 2012 """

    return 0.5 * abs(self.void_barrier) * (1/(2*pi)**0.5) \
           * (1/(sigma**3)) * exp((-1 * self.void_barrier**2)/(2*sigma**2))

  multiplicity_function_pls_vec = np.vectorize(multiplicity_function_pls)

  def scaled_void_distribution(self,nu):
    """Scaled distribution of void masses/sizes
       from Sheth & van de Weyagert (2004 MNRAS 350 517)

    Parameters
    ----------

    nu : array
    scaled void mass/size

    void_barrier : float
    the critical underdensity; defines a void

    collapse_barrier : float
    the critical overdensity; boundary at which virialized objects are formed

    Notes
    -----

    Calculates nu given in equation [4]; distribution demonstrated in figure(7)
    """

    # D ; the void-and-cloud parameter
    D = self.void_and_cloud(self.void_barrier,self.collapse_barrier)

    return self.multiplicity_function_jlh_vec(nu,D)


  def void_radii_dist(self,logr,cosm):
    """ Number density wrt void radii, Jennings, Li & Hu (2013 MNRAS 000 1)

    Parameters
    ----------

    r : array
    void radii

    ps : object
    PowSpec class instance

    f : function
    Multiplicity function

    z : float
    redshift

    Notes
    -----

    Produces the differential number density of voids
    wrt to their characteristic radius.
    equation [12], Jennings, Li & Hu
    """

    # account for expansion factor (eq.12, Jennings,Li & Hu)
    logr_l = log(exp(logr) / 1.7)

    # D ; the void-and-cloud parameter
    D = self.void_and_cloud()

    # calculate volume from a given R
    V = (4. * pi * pow(exp(logr),3)) / 3.

    # get sigma from PowSpec class fit
    sigma = cosm.sig_fit_z(logr)

    # get dln(1/sigma) / dln(r)
    dlns_dlnr = fabs(cosm.dlnsigmadlnr(logr_l))

    # calculate f(sigma)
    if(np.isscalar(sigma)):
      fSig = self.multiplicity_function_jlh(sigma,D)
    else:
      fSig = self.multiplicity_function_jlh_vec(self,sigma,D)

    no_dens = (fSig/V) * dlns_dlnr

    return no_dens

  def void_radii_dist_pls(self,r,cosm):
    """ Number density wrt void radii, Paranjape, Lam & Sheth 2011

    Parameters
    ----------

    r : array
    void radii

    ps : object
    PowSpec class instance

    f : function
    Multiplicity function

    z : float
    redshift

    Notes
    -----

    Produces the differential number density of voids
    wrt to their characteristic radius.
    Based on Jennings, Li & Hu (2013 MNRAS 000 1) equation [12]
    """

    # account for expansion factor (eq.12, Jennings,Li & Hu)
    r = (r / 1.7)

    # D ; the void-and-cloud parameter
    D = self.void_and_cloud()

    # calculate volume from a given R
    V = (4. * pi * pow(r,3)) / 3.

    # get sigma from PowSpec class fit
    sigma = cosm.Dplus * cosm.sig_fit(log(r*1.7))

    # get dln(1/sigma) / dln(r)
    dlns_dlnr = fabs(cosm.dlnsigmadlnr(log(r)))

    # calculate f(sigma)
    if(np.isscalar(sigma)):
      fSig = self.multiplicity_function_pls(sigma)
    else:
      fSig = self.multiplicity_function_pls_vec(self,sigma)

    no_dens = (fSig/V) * dlns_dlnr

    return no_dens

  def void_radii_dist_linear(self,r,cosm):
    """ Number density wrt void radii in the
        linear domain, Jennings, Li & Hu (2013 MNRAS 000 1)

    Parameters
    ----------

    r : array
    void radii

    cosm : object
    Cosmology class instance

    z : float
    redshift    pdf.append(evs.void_survey_pdf(i,norm,cosm,vd,10**5))



    Notes
    -----

    Produces the differential number density of voids
    wrt to their characteristic radius in the linear domain
    (equation [10] Jennings,Li & Hu)
    """

    # D ; the void-and-cloud parameter
    D = self.void_and_cloud()

    # calculate volume from a given R
    V = (4. * pi * pow(r,3)) / 3.

    # get sigma from PowSpec class fit
    sigma = cosm.sig_fit(log(r))

    # get dln(1/sigma) / dln(r)
    dlns_dlnr = fabs(cosm.dlnsigmadlnr(log(r)))

    # calculate f(sigma)
    if(np.isscalar(sigma)):
      fSig = self.multiplicity_function_jlh(self,sigma,D)
    else:
      fSig = self.multiplicity_function_jlh_vec(sigma,D)

    no_dens = (fSig/V) * dlns_dlnr

    return no_dens


  def void_radii_dist_vdn(self,r,cosm):
    """ Volume conserving model for void radii
    relative abundances, Jennings, Li & Hu 2013 """

    # D ; the void-and-cloud parameter
    D = self.void_and_cloud()

    # calculate volume from a given R
    V = (4. * pi * pow(r,3)) / 3.

    # get sigma from PowSpec class fit
    sigma = cosm.Dplus * cosm.sig_fit(log(r))

    # get dln(1/sigma) / dln(r)
    dlns_dlnr = fabs(cosm.dlnsigmadlnr(log(r/1.7)))

    # calculate f(sigma)
    if(np.isscalar(sigma)):
      fSig = self.multiplicity_function_jlh(sigma,D)
    else:
      fSig = self.multiplicity_function_jlh_vec(self,sigma,D)

    no_dens = (fSig/V) * dlns_dlnr

    return no_dens

  def void_radii_dist_vdn_pls(self,r,cosm):
    """ Volume conserving model for void radii
    relative abundances, Jennings, Li & Hu 2013

    Uses the correlated steps approach multiplicity
    function (Paranjape et al. 2012) """

    # D ; the void-and-cloud parameter
    D = self.void_and_cloud()

    # calculate volume from a given R
    V = (4. * pi * pow(r,3)) / 3.

    # get sigma from PowSpec class fit
    sigma = cosm.Dplus * cosm.sig_fit(log(r))

    # get dln(1/sigma) / dln(r)
    dlns_dlnr = fabs(cosm.dlnsigmadlnr(log(r/1.7)))

    # calculate f(sigma)
    if(np.isscalar(sigma)):
      fSig = self.multiplicity_function_pls(sigma)
    else:
      fSig = self.multiplicity_function_pls_vec(self,sigma)

    no_dens = (fSig/V) * dlns_dlnr

    return no_dens

  def void_radii_dist_ps(self,r,cosm):
    """ void abundance wrt radius from press schechter formalism
        (as defined in Kamionkowski, Verde & Jimenez (2008), equation [4])
    """
    rho = cosm.rho_m(0.0)

    m = (4 * pi * rho * r**3)/3

    sigma = cosm.Dplus * cosm.sig_fit(log(r))

    dlns_dlnm = fabs(cosm.dlnsigmadlnm(log(m)))

    n = (9/(2*pi**2)) * ((pi/2)**0.5) * (1/r**4) * (abs(self.void_barrier)/sigma) \
     * abs(dlns_dlnm) * exp(-1*(self.void_barrier**2/(2*sigma**2)))

    return n

  def void_radii_dist_ps_ng(self,r,S3,cosm):
    """ void abundance wrt radius from press schechter formalism with Non-Gaussian IC's
        (as defined in Kamionkowski, Verde & Jimenez (2008), equation [6])
    """

    rho = cosm.rho_m(0.0)

    m = (4 * pi * rho * r**3)/3

    # differential skewness wrt mass
    dS3dm = self.ds3_dm_fit(m)

    sigma = cosm.Dplus * cosm.sig_fit(log(r))

    dlns_dlnm = fabs(cosm.dlnsigmadlnm(log(m)))

    a = (9/(2*pi**2)) * ((pi/2)**0.5) * (1/r**4) * exp(-1*(self.void_barrier**2/(2*sigma**2)))
    b = (abs(self.void_barrier)/sigma) - ((S3*sigma)/6)*( ((self.void_barrier/sigma)**4) \
    - 2*((self.void_barrier/sigma)**2) - 1 )
    c = (1/6)*dS3dm*sigma*(((self.void_barrier/sigma)**2)-1)

    n = a * ( (dlns_dlnm*b) + c )

    return n

  def skewness_S3(self,r,cosm):
    """ Skewness parameter S3 from Enqvist, Hotchkiss & Taanila (2010)
        Equations [2.7] stated below; confusion as to the cosmology
        dependence of this form, since sigma_8 is stated
    """
    #P = cosm.power_spectrum_P(1/r,z)
    sigma = cosm.sig_fit(log(r))
    #s3 = (6 * cosm.fnl * (P**0.5)) / sigma
    s3 = (3. * 10**-4) * cosm.fnl / sigma

    return s3

  skewness_S3_vec = np.vectorize(skewness_S3)

  def ds3dm(self,m,s3):
    """ slope of skewness wrt mass scale
        dS3/dM
    """
    fit = np.poly1d(np.polyfit(m,s3,10))
    self.ds3_dm_fit = np.polyder(fit)

    return self.ds3_dm_fit

  def cumulative_V_R(self,r,cosm,func):
    """ cumulative fration of volume contained within voids
    from Jennings, Li & Hu, 2013 """

    def cumulative_int(r,cosm,func):
      # integral function argument
      n = func(r,cosm)
      V = (4 * pi * r**3) / 3
      return n * V * (1/r)

    return integrate.quad(cumulative_int,r,100.,args=(cosm,func))[0]

  cumulative_V_R_vec = np.vectorize(cumulative_V_R)

  def cumulative_V_R_NG(self,r,cosm,func):
    """ cumulative fration of volume contained within voids
    from Jennings, Li & Hu, 2013

    NON - GAUSSIAN """

    def cumulative_int(r,vd,cosm,func):
      s3 = vd.skewness_S3(r,cosm)
      # integral function argument
      n = func(r,s3,cosm)
      V = (4 * pi * r**3) / 3
      return n * V * (1/r)

    return integrate.quad(cumulative_int,r,np.inf,args=(self,cosm,func))[0]

  cumulative_V_R_NG_vec = np.vectorize(cumulative_V_R_NG)

  def void_mass_dist(self,m,cosm,z=0.0):
    """ Number density wrt void mass, Jennings, Li & Hu (2013 MNRAS 000 1)

    Parameters
    ----------

    m : array
    void masses

    cosm : object
    Cosmology class instance

    cosm : object
    Cosmology class instance

    z : float
    redshift

    Notes
    -----

    produces the differential number density of voids
    wrt their characteristic mass
    """
    # D ; the void-and-cloud parameter
    D = self.void_and_cloud()

    # Average density of the universe
    rho = cosm.rho_m(z)

    # convert mass -> radius,use to get sigma from PowSpec class fit
    r = ( ((3*m)/(4*pi*rho))**(0.3333333333333333333) ) / 1.7
    sigma = cosm.sig_fit(log(r*1.7))

    # get dln(1/sigma) / dln(m)
    dlns_dlnm = fabs(cosm.dlnsigmadlnm(log(m)))

    # calculate f(sigma)
    if(np.isscalar(sigma)):
      fSig = self.multiplicity_function_jlh(sigma,D)
    else:
      fSig = self.multiplicity_function_jlh_vec(sigma,D)

    no_dens = ((fSig*rho)/m) * dlns_dlnm

    #plt.loglog(m,fSig,m,dlns_dlnm)
    #plt.legend(["fSig","dlns_dlnm"])
    #plt.show()
    #raw_input()

    return no_dens

  def conditional_mf_bond(self,halo_m,void_m,cosm,z=0.0):
    """ Conditional Mass function, Furlanetto & Piran 2008
        equation [5]
        based on the Bond et al. prescripton of the Press Schechter mass function
    """
    rho = cosm.rho_m(z)

    r = ((halo_m*3)/(4*pi*rho))**(0.333333333333333)
    sigma = cosm.Dplus * cosm.sig_fit(log(r))

    # void radius for a region with linearized underdensity delta_v
    r = ((void_m*3)/(4*pi*rho*0.2))**(0.333333333333333)
    sigmav = cosm.Dplus * cosm.sig_fit(log(r))

    n =  ((2/pi)**0.5) * (rho/halo_m**2) * abs(cosm.dlnsigmadlnm(log(halo_m))) * (sigma**2) \
    * ((self.collapse_barrier-self.void_barrier)/((sigma**2) - (sigmav**2))**1.5) \
    * exp((-1 * (self.collapse_barrier-self.void_barrier)**2)/(2 * ((sigma**2) - (sigmav**2))))

    return n

  def massfunction_ps(self,halo_m,void_m,cosm,z=0.0):
    """ non-conditional mass function, based on Press & Schechter 1974
        equation [4] in Furlanetto & Piran
    """
    rho = cosm.rho_m(z)

    r = ((halo_m*3)/(4*pi*rho))**(0.333333333333333)
    sigma = cosm.Dplus * cosm.sig_fit(log(r))

    n = ((2/pi)**0.5) * (rho/(halo_m**2)) * abs(cosm.dlnsigmadlnm(log(halo_m))) \
        * (self.collapse_barrier/sigma) * exp((-self.collapse_barrier**2)/(2*(sigma**2)))

    return n

  def hod_kravtsov(self,halo_m,m_min):
    """ Halo occupation distribution (kravtsov et al. 2004)

    Parameters
    ----------

    halo_m : float
    Halo mass

    C : int
    Halo ocupation model fit parameter

    beta : int
    Halo occupation model fit parameter

    m_min : float
    minimum galaxy mass, as defined by the obervational detection threshold,
    NOT the actual minimum mass of a galaxy (Kravtsov et al. 2004)

    Notes
    -----

    As demonstrated in Furlanetto & Piran, equation [6]. Total occupation
    number composed of the central occupation plus the satellite occupation
    """
    # recommended values for z = 0 from Kravtsov et al.
    C = 30
    beta = 1

    # satellite galaxy occupation number
    Ns = (halo_m / (C * m_min))**beta

    # central galaxy occupation number
    Nc = 1.

    return Ns + Nc

  def galaxy_no_density(self,void_m,m_min,cosm,z=0.0,conditional=True):
    """ Total comoving number density of galaxies within a given region
        Furlanetto & Piran 2008, equation [7]
    """

    def integ(halo_m,void_m,m_min,cosm,z,conditional):
      if conditional is True:
        ng = self.hod_kravtsov(halo_m,m_min) * self.conditional_mf_bond(halo_m,void_m,cosm,z)
      else:
        ng = self.hod_kravtsov(halo_m,m_min) * self.massfunction_ps(halo_m,void_m,cosm,z)
      return ng

    n_gal, error = integrate.quad(integ,m_min,np.inf,args=(void_m,m_min,cosm,z,conditional))

    return n_gal

  def galaxy_underdensity(self,void_m,m_min,cosm,z=0.0):
    """ total observed galaxy underdensity in a void with physical size Rv
        Furlanetto & Piran 2008, equation [8]
    """
    n_gal = self.galaxy_no_density(void_m,m_min,cosm,z)
    n_gal_ps = self.galaxy_no_density(void_m,m_min,cosm,z,conditional=False)

    delta_g = (n_gal/(n_gal_ps * (1.7**3)))-1
    return delta_g

  def galaxy_ud_fit(self,delta_g,delta_v):
    """ Calculates a fit to the galaxy underdensity / void underdensity relationship
    """
    # fit a polynomial to the given data
    fit = np.polyfit(delta_v,delta_g,4)
    # differentiate and solve to find maximum along delta_v
    dv = np.polyder(np.poly1d(fit))
    max_fit = np.roots(dv)

    print "max_fit: %s" % max_fit

    # find the maximum within the specified range
    for x in max_fit:
      if min(delta_v) < x < max(delta_v):
        max_fit = x
        max_index = (np.abs(delta_v-max_fit)).argmin()
        break

    dv_new = []
    dg_new = []
    #reduce original arrays between 0 and maximum
    for i,x in enumerate(delta_v):
      if i < max_index: dv_new.append(x)
    for i,x in enumerate(delta_g):
      if i < max_index: dg_new.append(x)

    #fit to function within range (avoids multiple solutions)
    self.delta_g_fit = np.poly1d(np.polyfit(dg_new,dv_new,10))

    return self.delta_g_fit

  def halo_underdensity(self,halo_m,void_m,cosm,z=0.0):
    """ Halo underdensity at mass m
    """
    n_halo = self.conditional_mf_bond(halo_m,void_m,cosm,z)
    n_halo_ps = self.massfunction_ps(halo_m,void_m,cosm,z)

    delta_h = (n_halo/(n_halo_ps*(1.7**3)))-1
    return delta_h

  def halo_ud_fit(self,delta_h,delta_v):
    """ Calculates a fit to the halo underdensity / void underdensity relationship
    """
    fit = np.polyfit(delta_h,delta_v,6)
    dv = np.poly1d(fit)

    self.delta_h_fit = np.poly1d(np.roots(fit))
    print self.delta_h_fit
    return self.delta_h_fit


def initialisation(radius,mass,cosm,z=0.0):
  cosm.vd_initialisation(radius,mass,z)
  #m_radius = ( (3 * mass) / (4 * pi * cosm.rho_m(z)) )**(0.33333333333)
  #vd.ds3dm(mass,skewness_S3(m_radius,cosm,z))

  dump_pickle(cosm)
  return True

def dump_pickle(cosm):
  f = open("void_init.p","wb")
  pickle.dump([cosm.Dplus,cosm.sigmar,cosm.sigmar_z,cosm.sig_fit,cosm.sig_fit_z,cosm.dlnsigmadlnr,cosm.dlnsigmadlnm],f)
  f.close()
  return True

def load_pickle():
  f = open("void_init.p","rb")
  Dz, sigmar, sigmar_z, sig_fit, sig_fit_z, dlnsigmadlnr, dlnsigmadlnm= pickle.load(f)
  f.close()
  return Dz, sigmar, sigmar_z, sig_fit, sig_fit_z, dlnsigmadlnr, dlnsigmadlnm

if __name__ == '__main__':
  cosm = Cosmology()
  vd = Void(cosm)

  radius = np.logspace(-1,2,600)
  mass = np.logspace(8.1,17,700)
  radius.tolist()
  mass.tolist()

  #initialisation(radius,mass,cosm,z=0.)

  # load pickled void parameters from powspec class
  cosm.Dplus, cosm.sigmar, cosm.sigmar_z, cosm.sig_fit, cosm.sig_fit_z, cosm.dlnsigmadlnr, cosm.dlnsigmadlnm = load_pickle()

  """
  #Implementation of Non-Gaussian methods
  # set non-gaussianity as non-zero
  cosm.fnl = -200.

  m_radius = ( (3. * mass) / (4. * pi * cosm.rho_m(0.)) )**(0.33333333333)
  r_mass = (4*pi*cosm.rho_m(0.)*radius**3)/3.

  # calculate the skewness for given set of radii
  s3 = []
  for r in radius:
    s3.append(vd.skewness_S3(r,cosm,cosm.pk))

  #plt.loglog(radius,s3)
  #plt.show()
  #raw_input()

  vd.ds3dm(r_mass,s3)

  n = vd.void_radii_dist_ps(radius,cosm.pk)
  n_NG = vd.void_radii_dist_ps_ng(radius,s3,cosm.pk)

  plt.title("f_nl = -200")
  plt.loglog(radius,n,radius,n_NG)
  plt.legend(["Gaussian","NG"])
  plt.xlabel("R")
  plt.ylabel("n")
  plt.show()

  raw_input()

  plt.semilogx(radius,n_NG/n,)
  plt.xlabel("R")
  plt.ylabel("ratio")
  plt.show()
  """

  """
  M_void = 10**14
  m_min = 2.5 * 10**11

  delta_g = []
  delta_h = []
  delta_v = np.arange(-6.0,0.0,0.01)
  for dv in delta_v:
    vd.void_barrier = dv
    #delta_h.append(vd.halo_underdensity(10**11,10**14,cosm,cosm.pk))
    delta_g.append(vd.galaxy_underdensity(M_void,m_min,cosm,cosm.pk))

  r = ((3*M_void)/(4*pi*cosm.rho_m(0.0)))**(0.3333333333333)
  sigma_2 = (cosm.pk.sig_fit(log(r)))**2

  vd.galaxy_ud_fit(delta_g,delta_v)
  #vd.halo_ud_fit(delta_h,delta_v)

  print "delta_v : %s || sigma_2 : %s" % (vd.delta_g_fit(-0.8), sigma_2)

  #plt.plot(delta_h,vd.delta_h_fit(delta_h),delta_h,delta_v)
  plt.plot(delta_g,vd.delta_g_fit(delta_g),delta_g,delta_v)
  #plt.plot(dg_new,vd.delta_g_fit(dg_new),dg_new,dv_new)

  plt.legend(["fit","real"])
  plt.ylabel(r"$\delta^{L}_{v}$")
  plt.xlabel("$\delta_{g}$")
  plt.show()

  vd.void_barrier = -2.7

  """

  """
  no_dens = []
  vd.delta_g_fit(-0.8)
  no_dens.append(vd.void_radii_dist(radius,cosm.pk))


  plt.loglog(radius,no_dens)
  plt.show()
  """

  """
  no_mass = []
  for m in mass:
    #no_mass.append(
    #print (vd.conditional_mf_bond(m,10**15,cosm,cosm.pk))
    #vd.void_barrier = 0.00
    #no_mass2.append(vd.conditional_mf_bond(m,0.,cosm,cosm.pk,sigmav=0.))
    #no_mass.append(
    #print (vd.hod_kravtsov(m,10**10))
    V = m / cosm.rho_m(0.0)
    no_mass.append(V * vd.galaxy_no_density(m,10**11,cosm,cosm.pk))

  plt.loglog(mass,no_mass)
  plt.title("Cumulative Galaxy no against resident Void Mass")
  plt.xlabel(r"Galaxy No.")
  plt.ylabel(r"Void Mass M_{0}")
  #plt.yscale((10**-12),1.)
  plt.show()
  """

  """
  halo_mass = np.logspace(7,17.5,700)
  halo_mass.tolist()

  gal_nd = []
  gal_nd2 = []
  gal_nd3 = []
  gal_nd4 = []

  ax1 = fig.add_subplot(111)

  for m in halo_mass:
    gal_nd.append(vd.conditional_mf_bond(m,10**16,cosm,cosm.pk))
    gal_nd2.append(vd.conditional_mf_bond(m,10**15,cosm,cosm.pk))
    gal_nd3.append(vd.conditional_mf_bond(m,10**14,cosm,cosm.pk))
    gal_nd4.append(vd.conditional_mf_bond(m,10**13,cosm,cosm.pk))

  #no_mass = vd.galaxy_no_density(10**16,10**10,cosm,cosm.pk)

  #print no_mass
  ax1.loglog(halo_mass,gal_nd)
  ax1.loglog(halo_mass,gal_nd2)
  ax1.loglog(halo_mass,gal_nd3)
  ax1.loglog(halo_mass,gal_nd4)

  ax1.legend([r"10^{16}",r"10^{15}",r"10^{14}",r"10^{13}"])

  plt.title("Comoving No Density of Galaxies against resident Halo Mass, \n in a Void of given Mass")
  plt.xlabel(r"Halo mass M_{0}")
  plt.ylabel("Galaxy No. Density")
  #plt.yscale((10**-12),1.)
  plt.show()
  """

  """
  #plot
  cumul, error = vd.cumulative_V_R_vec(vd,radius,cosm.pk,vd.void_radii_dist_vdn)
  vd.collapse_barrier = 1.686
  cumul2, error = vd.cumulative_V_R_vec(vd,radius,cosm.pk,vd.void_radii_dist_vdn)
  plt.fill_between(radius,cumul,cumul2,alpha=0.5)

  vd.collapse_barrier = 1.06
  cumulb, error = vd.cumulative_V_R_vec(vd,radius,cosm.pk,vd.void_radii_dist_linear)
  vd.collapse_barrier = 1.686
  cumulb2, error = vd.cumulative_V_R_vec(vd,radius,cosm.pk,vd.void_radii_dist_linear)
  plt.fill_between(radius,cumulb,cumulb2,alpha=0.5)

  vd.collapse_barrier = 1.06
  cumulc, error = vd.cumulative_V_R_vec(vd,radius,cosm.pk,vd.void_radii_dist)
  vd.collapse_barrier = 1.686
  cumulc2, error = vd.cumulative_V_R_vec(vd,radius,cosm.pk,vd.void_radii_dist)
  plt.fill_between(radius,cumulc,cumulc2,alpha=0.5)

  plt.xscale('log')
  plt.show()

  """

  """
  r_mass = (4*pi*cosm.rho_m(0.)*radius**3)/3.

  # calculate the skewness for given set of radii
  s3 = []
  cosm.fnl = 200.
  s3 = (vd.skewness_S3_vec(vd,radius,cosm))

  # differential swekness wrt mass
  vd.ds3dm(r_mass,s3)

  #ps_ng_nr = vd.void_radii_dist_ps_ng(radius,s3,cosm)
  #cumul_NG = vd.cumulative_V_R_NG_vec(vd,radius,cosm,vd.void_radii_dist_ps_ng)
  cosm.fnl = 0.


  # calculate the number density of voids wrt radius
  #no_rad = vd.void_radii_dist(radius,cosm)
  #no_rad_pls = vd.void_radii_dist_pls(radius,cosm)
  #no_radius = vd.void_radii_dist_vdn(radius,cosm)
  #no_radius2 = vd.void_radii_dist_vdn_pls(radius,cosm)
  #ps_nr = vd.void_radii_dist_ps(radius,cosm)


  # caculate the cumulative volume fraction in voids with radius r
  cumul = vd.cumulative_V_R_vec(vd,radius,cosm,vd.void_radii_dist)
  cumul2 = vd.cumulative_V_R_vec(vd,radius,cosm,vd.void_radii_dist_pls)
  cumul3 = vd.cumulative_V_R_vec(vd,radius,cosm,vd.void_radii_dist_vdn)
  cumul4 = vd.cumulative_V_R_vec(vd,radius,cosm,vd.void_radii_dist_vdn_pls)
  cumul5 = vd.cumulative_V_R_vec(vd,radius,cosm,vd.void_radii_dist_ps)

  #plt.loglog(radius,no_rad,radius,no_rad_pls,radius,no_radius,radius,no_radius2,radius,ps_nr,radius,ps_ng_nr)
  plt.semilogx(radius,cumul,radius,cumul2,radius,cumul3,radius,cumul4,radius,cumul5)
  plt.legend(["SVDW","PLS","JLH","JLH-PLS","PS"])

  #plt.xlim([1e-1,5e1])
  #plt.ylim([1e-10,1e-0])

  plt.ylim(0.,2.)
  plt.xlim(1e-1,5e1)


  plt.xlabel(r'r [Mpc/h]', fontsize='15')
  #plt.ylabel(r'dn/dlnr $(h/Mpc)^{3}$', fontsize='15')
  plt.ylabel(r"F(R)",fontsize='15')

  plt.show()

  """


  """
  norm = evs.void_norm(vd,cosm,vd.void_radii_dist)
  pdf = evs.void_pdf_vec(radius,norm,cosm,vd,vd.void_radii_dist,10**6)

  norm = evs.void_norm(vd,cosm,vd.void_radii_dist_pls)
  pdf2 = evs.void_pdf_vec(radius,norm,cosm,vd,vd.void_radii_dist_pls,10**6)

  norm = evs.void_norm(vd,cosm,vd.void_radii_dist_vdn)
  pdf2 = evs.void_pdf_vec(radius,norm,cosm,vd,vd.void_radii_dist_vdn,10**6)

  norm = evs.void_norm(vd,cosm,vd.void_radii_dist_vdn_pls)
  pdf2 = evs.void_pdf_vec(radius,norm,cosm,vd,vd.void_radii_dist_vdn_pls,10**6)

  norm = evs.void_norm(vd,cosm,vd.void_radii_dist_ps)
  pdf2 = evs.void_pdf_vec(radius,norm,cosm,vd,vd.void_radii_dist_ps,10**6)

  plt.xlim(1e1,5e1)
  plt.ylim(1e-7,1e-0)
  plt.loglog(radius,pdf,raius,pdf2)
  plt.show()

  """

  """

  nod = []
  nod2 = []

  nu_range = np.arange(0.1,20,0.1)

  for R in nu_range:
    nod.append(vd.void_radii_dist(R,void_barrier=-2.7,cosm.pk))
    nod2.append(vd.void_radii_dist(R,collapse_barrier=1.686,void_barrier=-2.7,cosm.pk))

  plt.plot(nu_range,nod,nu_range,nod2)

  plt.yscale('log')
  plt.xscale('log')
  #plt.yscale('linear')
  #plt.xscale('linear')

  plt.xlabel(r'r [Mpc/h]', fontsize='20')
  plt.ylabel(r'dn/dlnr $(h/Mpc)^{3}$', fontsize='20')
  plt.legend((r'$\delta_{c}=1.06$',r'$\delta_{c}=1.69$'), prop={'size':20})

  plt.show()
  """
  """
  nu_range = np.arange(0.1,20,0.1)

  fnu1 = []
  fnu2 = []
  fnu3 = []

  D = void_and_cloud(void_barrier=-2.81,collapse_barrier=1.06)

  for nu in nu_range:
    fnu1.append(scaled_void_distribution(nu,collapse_barrier=1.06))
    fnu2.append(scaled_void_distribution(nu,collapse_barrier=1.69))
    fnu3.append(scaled_void_distribution(nu,collapse_barrier=99999999))

  plt.plot(nu_range,fnu1,nu_range,fnu2,nu_range,fnu3)
  plt.legend((r'$\delta_{c}=1.06$',r'$\delta_{c}=1.69$',r'$\delta_{c}=\infty$'), prop={'size':20})
  plt.ylim(0.01,0.6)
  plt.xlim(0.,5)
  #plt.xscale('log')
  #plt.yscale('log')
  plt.yscale('linear')
  plt.xscale('linear')
  plt.xlabel(r'$\nu = (\delta_{V} / \sigma)^{2}$', fontsize='20')
  plt.ylabel(r'$f(\nu)$', fontsize='20')

  plt.show()
"""
