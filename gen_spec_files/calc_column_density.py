# Calculates the column density for a molecule using fixed mass mixing ratio
# from the top of the atmosphere down to a P_surf pressure surface.

from pylab import *
from netCDF4 import *

# Add the opacity_tools/cmn folder to the Python path
import sys
import os
sys.path.append(os.path.abspath('../cmn'))

# Import cmn modules
from gas_list import gas_name

# Do not expect ATMO abundance file by default
use_atmo_chem_file = False

################################################################################
# Planet parameters
################################################################################

# HD 209458b
# from mass_mixing_ratio_solar import mass_mixing_ratio, max_mass_mixing_ratio
# from mass_mixing_ratio_atmo import mass_mixing_ratio, use_atmo_chem_file
# grav_acc = 13.9 # Acceleration of gravity
# P_surf = 1e+8 # Column density is calculated down to P_surf
# P_ref = 1e+4 # Reference pressure at which to use abundance
# T_ref = 1e+3 # Reference temperature at which to use abundance
# flag = ''
# gases = array([
#     'H2',
#     'He',
#     'CO',
#     'CH4',
#     'H2O',
#     'NH3',
#     'TiO',
#     'VO',
#     'Na',
#     'K',
#     'Li',
#     'Rb',
#     'Cs',
#     'FeH',
#     'CrH',
#     'PH3',
#     'C2H2',
#     'HCN',
#     'H2S',
#     ])
# atmo_chem_file = 'chem_hd189_eq.ncdf'

# Archean Earth
from mass_mixing_ratio_archean import (mass_mixing_ratio,
    max_mass_mixing_ratio, p_std)
grav_acc = 9.81 # Acceleration of gravity
P_surf = 1e+5 # Column density is calculated down to P_surf
P_ref = 1e+4 # Reference pressure at which to use abundance
T_ref = 273.15 # Reference temperature at which to use abundance
flag = 'Charnay_et_al_Case_A' # Flag provided to mixing ratio routine
gases = array([
    'H2O',
    'CO2',
    'CH4',
    ])
integrate_column = True

################################################################################
# End user input
################################################################################

# Print flag if provided
if flag:
  print('Using mass mixing ratio flag ' + flag)
  print()

# Calculate maximum column density
if not use_atmo_chem_file:
  print('Maximum column masses')
  for gas in gases:
    if flag:
      column_density = max_mass_mixing_ratio(gas,flag=flag)*P_surf/grav_acc
    else:
      column_density = max_mass_mixing_ratio(gas)*P_surf/grav_acc
    print(gas.ljust(10), '{0:1.1e}'.format(column_density))

# Calculate column density using abundances at P_ref, T_ref
print()
print('Column masses using abundances at ' +
    'P_ref = {0:1.1e} Pa, T_ref = {1:g} K'.format(P_ref, T_ref))
for gas in gases:
  if use_atmo_chem_file:
    column_density = (mass_mixing_ratio(atmo_chem_file, gas, P_ref)*
        P_surf/grav_acc)
  else:
    if flag:
      column_density = mass_mixing_ratio(gas, P_ref, T_ref,
          flag=flag)*P_surf/grav_acc
    else:
      column_density = mass_mixing_ratio(gas, P_ref, T_ref)*P_surf/grav_acc
  print(gas.ljust(10), '{0:1.1e}'.format(column_density))

print()
print('Vertically integrated column masses')
for gas in gases:
  if use_atmo_chem_file:
    column_mass = sum(mass_mixing_ratio(gas, (p_std[1:] + p_std[:-1])/2.0)*
          diff(-p_std)/grav_acc)
  else:
    if flag:
      column_mass = sum(mass_mixing_ratio(gas, (p_std[1:] + p_std[:-1])/2.0,
          T_ref, flag=flag)*diff(-p_std)/grav_acc)
    else:
      column_mass = sum(mass_mixing_ratio(gas, (p_std[1:] + p_std[:-1])/2.0,
          T_ref)*diff(-p_std)/grav_acc)
  print(gas.ljust(10), '{0:1.1e}'.format(column_mass))