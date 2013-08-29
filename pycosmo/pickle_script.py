"""
Pickle script, for creating, finding and loading pickles in PyCosmo
"""

import os
from os import listdir
from os.path import isfile, join

import cPickle as pickle
import utils as ut

import cosmology

def init_pickle(cosm,default):
  pick = search_pickle()
  print "Choose a previously defined cosmology, \n or enter '0' to import a new cosmology: \n"

  for i,x in enumerate(pick):
    print "[%s] - %s" % (i+1,x)

  while True:
    try:
      a = int(raw_input("\n:"))
      break
    except ValueError:
      print "Invalid entry, please try again"

  if a is 0:
    initialise_cosmology(cosm,default)
  else:
    load_pickle(str(pick[a-1]),cosm)

  return True

def search_pickle():
  """ search for pickles available """
  mypath = os.path.join(os.path.dirname(os.path.realpath(__file__)),"pickle/")
  files = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]
  return files

def load_pickle(name,cosm):
  """ load current pickle """
  mypath = os.path.join(os.path.dirname(os.path.realpath(__file__)),"pickle/")
  f = open(os.path.join(mypath,name),"rb")

  cosm.delta_c,cosm.h_0,cosm.O_m0,cosm.O_de0,cosm.w_0,cosm.O_b0,cosm.O_r0,cosm.O_k0,cosm.fnl,cosm.tau_rec,cosm.z_rec,cosm.H_0,cosm.n_s=pickle.load(f)

  #a = pickle.load(f)
  #print a

  #cosm.Dz, cosm.sigmar, cosm.sig_fit, cosm.dlnsigmadlnr, cosm.dlnsigmadlnm = pickle.load(f)
  f.close()
  return True

def dump_pickle(cosm,name,default):
  """ create new pickle """

  mypath = os.path.join(os.path.dirname(os.path.realpath(__file__)),"pickle/")
  f = open(os.path.join(mypath,name),"wb")

  pickle.dump([default["dc"],default["h0"],default["om"],default["ode"],default["w0"],default["ob"],default["o_r"],default["o_k"],default["f_nl"],default["tau_r"],default["z_r"],default["h0"]*100.e0,default["ns"]],f)

  """
  pickle.dump(default["h0"],f)
  pickle.dump(default["om"],f)
  pickle.dump(default["ode"],f)
  pickle.dump(default["w0"],f)
  pickle.dump(default["ob"],f)
  pickle.dump(default["o_r"],f)
  pickle.dump(default["o_k"],f)
  pickle.dump(default["f_nl"],f)
  pickle.dump(default["tau_r"],f)
  pickle.dump(default["z_r"],f)
  pickle.dump(default["h0"]*100.e0,f)
  pickle.dump(default["ns"],f)
  """

  #pickle.dump([cosm.Dz,cosm.sigmar,cosm.sig_fit,cosm.dlnsigmadlnr,cosm.dlnsigmadlnm],f)
  f.close()

  load_pickle(name,cosm)

  return True

def initialise_cosmology(cosm,default):
  print "The following default parameters for a \n \
  WMAP7+BAO+H0 cosmology have been defined: \n"

  print default

  print "\n To change these parameters, open \n \
  and edit the 'default' dictionary within \n \
  cosmology.py"

  raw_input("If you are happy with these parameters \n \
  a pickle with further dependent variables \n \
  will now be created: [Enter to continue] \n \n ")

  while True:
    try:
      name = str(raw_input("Enter name of new cosmology:"))
      name += ".p"
      break
    except ValueError:
      print "Invalid entry, please try again"

  dump_pickle(cosm,name,default)

  return None
























