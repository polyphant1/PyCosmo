"""
Power spectrum script

Creates a power spectrum for given k
from either a CAMB matter power spectrum
or an analytic function calculated using
Eisenstein & Hu's fitting function
"""


import numpy as np
from numpy import sqrt, log, exp, fabs, pi, cos, sin

from scipy import integrate
from scipy import interpolate

import EH.power as power
import constants as ct


def powspec_initialise(cosm):
  if choose():
    z = raw_input("Enter z value of produced spectrum:")
    ident = raw_input("Enter identification number:")
    k,P = import_powerspectrum(ident,z)
    return interpolate(k,P)

def choose():
  """ Choose between an Eisenstein & HU fitting function or a CAMB power spectrum """
  print "  \n \
    ------------------------------------------------------------- \n \
    Would you like to import a CAMB matter/power spectrum file? \n \
    If not, an analytic spectrum using the prescription of \n \
    Eisenstein & Hu (1999,511) will be used. \n \n \
    If using a CAMB file ensure the cosmological \n \
    parameters used to generate the spectrum are \n \
    identical to those specified within the \n \
    cosmology class. \n \n \
    ------------------------------------------------------------- \n \
    If CAMB, enter 'True', else EH \n \
    ------------------------------------------------------------- \n"

  s = raw_input(":")

  if s == "True":
    return True
  return False

def import_powerspectrum(ident,z=0.0):
  """ import power spectrum function from
      a CAMB produced output file """

  z=str(int(z))

  k_array = []
  P_array = []

  for line in open('camb/camb_{0}_matterpower_z{1}.dat'.format(ident,z),'r'):
    a, b = line.split()
    k_array.append(float(a))
    P_array.append(float(b))

  return k_array, P_array

def interpolate(array_1,array_2):
  """ returns a function that uses interpolation
      to find the value of new points """

  return interpolate.interp1d(array_1,array_2)


def transfer_function_EH(cosm,k,z):
  """Calculates transfer function given wavenumber"""

  # set cosmology
  power.TFmdm_set_cosm(cosm.O_m0,cosm.O_b0,\
                      0.0,0,cosm.O_de0,cosm.h_0,z)

  """Call to EH power.c script
     ???? h Mpc^-1 OR Mpc^-1 ???? """
  #return power.TFmdm_onek_mpc(k)

  return power.TFmdm_onek_hmpc(k)


def power_spectrum_P(cosm,k,z):
  """ returns the power spectrum P(k) for a given z"""

  delta_h = 1.94 * (10**-5) * cosm.O_m0**(-0.785 - (0.05*log(cosm.O_m0))) \
             * exp(-0.95*(cosm.n_s-1)-0.169*(cosm.n_s-1)**2)

  Tk = transfer_function_EH(cosm,k,z)

  c_l = ct.const["c"] / ct.convert["Mpc_m"]   # speed of light in Mpc s^-1

  return (delta_h**2 * 2. * pi**2. * k**cosm.n_s) * \
         (c_l/(cosm.h_0 * ct.convert['H0']))**(3.+cosm.n_s) * (Tk*cosm.Dplus)**2
