# Functions to get P-T profiles and heating/cooling rate profiles from
# WRITER output in ModelE.

from pylab import *
import os

# Logical determining if the pt_profile routine defines the pressure
defines_pressure = True


# Reads a modelE PT profile from file, converts it to SI units, derives both
# level and layer values, and returns it
def pt_profile(profile_name):
  
  # Read file with data
  data = loadtxt(profile_name, skiprows=2, usecols=(1,3))
  
  # Extract level pressures (convert for mbar to Pa) and temperatures
  p_lev = data[:,0]*1e+2
  t_lev = data[:,1]
  
  # Calculate layer pressure and temperatures
  p_lay = (p_lev[:-1] + p_lev[1:])/2.0
  t_lay = (t_lev[:-1] + t_lev[1:])/2.0
  
  # Return everything
  return p_lev, t_lev, p_lay, t_lay


# Reads the (upwards and downwards) thermal flux from modelE and returns it
def read_thermal_flux(profile_name):

  # Read file with data
  data = loadtxt(profile_name, skiprows=2, usecols=(1,5,6))
  
  # Extract level pressures (convert for mbar to Pa) and fluxes
  p_lev = data[:,0]*1e+2
  flux_down = data[:,1]
  flux_up = data[:,2]
  
  # Return everything
  return p_lev, flux_down, flux_up
