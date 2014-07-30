"""
SDSS DR7 Void catalogue parsing script
"""

import numpy as np
from matplotlib import pyplot as plt
from numpy import sqrt, log, exp, fabs, pi, sin, cos

survey_volume_norm = []
survey_radius = []
survey_volume = []

#open and parse catalogue file
with open("../catalogues/centers_central_lss.dr72bright1.dat.out","r") as data:
  #skip first line
  data.readline()
  #parse rest of file
  for line in data.readlines():
    line = line.split()
    survey_volume_norm.append(float(line[3]))  # Volume (normalised)
    survey_radius.append(float(line[4]))       # Radius (Mpc/h)
    survey_volume.append(float(line[6]))       # Volume (Mpc/h^3)

# logarithmic bin space
bins = np.logspace(-0.5,1.2,30)

norm_radius = []

for x in survey_volume_norm:
  norm_radius.append( ((x * 3.)/ (4.*pi))**(0.333333333) )



max_survey_r = max(norm_radius)

"""

# bin survey data using numpy.histogram
freq, bins2 = np.histogram(norm_radius,bins)

# append extra value to account for smaller frequency array (allows plotting)
freq = np.append(freq,0)

# number density in (h/Mpc)**3
freq = freq / (1.9 * 10**6)

print freq, bins2

plt.loglog(bins2,freq)
plt.show()

"""
