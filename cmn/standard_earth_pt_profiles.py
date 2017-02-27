# Functions to get P-T profiles in standard Earth atmospheres.

from pylab import *
from scipy.interpolate import interp1d
import sys

# Logical determining if the pt_profile routine defines the pressure
defines_pressure = True

# Returns the selected P-T profile
def pt_profile(profile_name):

  data = loadtxt(profile_name, skiprows = 1, usecols = (1, 2))
  
  p_lev = data[:,0]*1.0e+2
  t_lev = data[:,1]
  
  # Pressures and temperatures in layers
  p_lay = (p_lev[:-1] + p_lev[1:])/2.0
  t_lay = (t_lev[:-1] + t_lev[1:])/2.0
  
  # Return everything
  return p_lev, t_lev, p_lay, t_lay