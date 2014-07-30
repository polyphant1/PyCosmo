"""
Cosmological EVS
"""

import numpy as np
from numpy import sqrt, log, exp, fabs, pi, sin, cos
from scipy import integrate

import void_distribution as vd


def void_fr(norm,r,cosm,vd,func):
  """ Original void distribution probability density function (pdf)

  Parameters
  ----------

  norm : float
  normalisation factor

  r : array
  void radii

  ps : object
  PowSpec class instance

  Notes
  -----

  f(r) from Harrison & Coles (2012); pdf of the original void distribution
  """

  dndlogr = func(r,cosm)
  return (1/norm)*dndlogr*(1/r)

def void_Fr(norm,r,cosm,vd,func,max_record):
  """ Cumulative void distribution

  Parameters
  ----------

  norm : float
  normalisation factor

  r : array
  void radii

  ps : object
  PowSpec class instance

  max_record : bool
  direction of cumulative sum

  Notes
  -----

  F(r) from Harrison & Coles (2012); known distribution of void radii

  for max_record=True, calculates the cumulative  distribution *upto*
  a given radii, otherwise calculates EVS for small radii
  """

  if max_record:
    integ = integrate.quad(func,0.,r,args=(cosm))[0]
  else:
    integ = integrate.quad(func,r,np.inf,args=(cosm))[0]

  return (1/norm)*integ

def void_norm(vd,cosm,func):
  """ Void Normalisation factor

  Parameters
  ----------

  ps : object
  PowSpec class instance

  Notes
  -----

  n_tot from Harrison & Coles (2012)
  normalisation factor; gives the
  total comoving number density of voids """

  return integrate.quad(func,0.,np.inf,args=(cosm))[0]

def void_pdf(r,norm,cosm,vd,func,V,max_record=True):
  """ Void probability density function (pdf)

  Parameters
  ----------

  r : array
  void radii

  norm : float
  normalisation factor

  ps : object
  PowSpec class instance

  V : float
  Constant redshift box volume

  max_record : bool
  direction of cumulative sum

  Notes
  -----

  phi(max) from Harrison & Coles (2012); exact extreme value pdf of the
  original void distribution for a given radius

  max_record passed to void_Fr function. Further information on its
  use provided there

  void pdf is redshift independent, i.e. also used for survey distributions
  """

  fr = void_fr(norm,r,cosm,vd,func)
  Fr = void_Fr(norm,r,cosm,vd,func,max_record)
  N = norm * V

  return N * fr * Fr**(N-1)

void_pdf_vec = np.vectorize(void_pdf)


def void_survey_fr(logr,norm,cosm,func,fsky,z_min,z_max):
  """ Original void distribution probability density function (pdf)
  in a given cosmological survey between min and max redshift values

  Parameters
  ----------

  vd : object
  Void class instance

  norm : float
  normalisation factor

  fsky : float
  sky fraction

  r : array
  void radii

  Notes
  -----

  f(r) from Harrison & Coles (2012); pdf of the original void distribution
  for a given survey (between mass and redshift limits)
  """

  def dndlnrdvdz(z,logr,cosm,func):
    cosm.growth(z)
    return func(logr,cosm) * cosm.dvdz(z)

  integ = integrate.quad(dndlnrdvdz,z_min,z_max,args=(logr,cosm,func))[0]

  return (fsky/norm)*integ

def void_survey_Fr(logr,norm,cosm,func,fsky,z_min,z_max,log_r_min):

  # integral function
  def dndlnrdvdz(rad,z,cosm,func):
    cosm.growth(z)
    return func(rad,cosm) * cosm.dvdz(z)

  #logr = np.linspace(log(r_min),log(r_max),r_steps)

  integ = integrate.dblquad(dndlnrdvdz, z_min, z_max, lambda rad: log_r_min, lambda rad: logr, args=(cosm,func))[0]

  # return F(r) values along with corresponding logr
  return (fsky/norm)*integ


def void_survey_norm(cosm,func,z_min=0.0,z_max=0.05,volume=1.9*(10**6),log_r_min=log(0.1),log_r_max=log(200)):
  fsky = volume / cosm.V_between(z_min,z_max)
  print fsky

  #integral function
  def dndlnrdvdz(logr,z,cosm,func):
     cosm.growth(z)
     return func(logr,cosm) * cosm.dvdz(z)

  integ = integrate.dblquad(dndlnrdvdz, z_min, z_max,
                            lambda r: log_r_min, lambda r: log_r_max, args=(cosm,func))[0]

  return fsky*integ

def void_survey_pdf(r,N,cosm,func,z_min=0.0,z_max=0.05,volume=1.9*(10**6),log_r_min=log(0.1)):
  fsky = volume / cosm.V_between(z_min,z_max)
  print fsky

  fr = void_survey_fr(log(r),N,cosm,func,fsky,z_min,z_max)
  Fr = void_survey_Fr(log(r),N,cosm,func,fsky,z_min,z_max,log_r_min)

  phi = N * fr * Fr**(int(N-1.))
  #logphi = log(N) + log(fr) + ((N-1.)*log(Fr))

  print r, fr, Fr, N, phi
  #print log(N), log(fr), ((N-1.)*log(Fr))

  return fr, Fr, phi

void_survey_pdf_vec = np.vectorize(void_survey_pdf)

