# Contains a function to calculate mass mixing ratios for various cases suitable
# for Paleo Mars.

from pylab import *
from scipy.interpolate import interp1d
from sys import exit

################################################################################
# Data
################################################################################

# External files with data required
from constants import r_gas, r_gas_cal
from gas_list import molar_weight

# Mean molecular weight of atmosphere set by call to mass_mixing_ratio
mean_molar_weight = None

# List of possible flags
mean_molar_weight = {
    'setup_1': 28.04,
    'setup_2':  2.56,
    'setup_3': 15.31,
    'setup_4': 39.53,
    }

# U.S. standard atmosphere 1976
# p_std = array([
#     1.013e+03, 8.988e+02, 7.950e+02, 7.012e+02, 6.166e+02,
#     5.405e+02, 4.722e+02, 4.111e+02, 3.565e+02, 3.080e+02,
#     2.650e+02, 2.270e+02, 1.940e+02, 1.658e+02, 1.417e+02,
#     1.211e+02, 1.035e+02, 8.850e+01, 7.565e+01, 6.467e+01,
#     5.529e+01, 4.729e+01, 4.047e+01, 3.467e+01, 2.972e+01,
#     2.549e+01, 1.743e+01, 1.197e+01, 8.258e+00, 5.746e+00,
#     4.041e+00, 2.871e+00, 2.060e+00, 1.491e+00, 1.090e+00,
#     7.978e-01, 4.250e-01, 2.190e-01, 1.090e-01, 5.220e-02,
#     2.400e-02, 1.050e-02, 4.460e-03, 1.840e-03, 7.600e-04,
#     3.200e-04, 1.450e-04, 7.100e-05, 4.010e-05, 2.540e-05,
#     ])*1.0e+2
# h2o_vmr_std = array([
#     7.745e+03, 6.071e+03, 4.631e+03, 3.182e+03, 2.158e+03,
#     1.397e+03, 9.254e+02, 5.720e+02, 3.667e+02, 1.583e+02,
#     6.996e+01, 3.613e+01, 1.906e+01, 1.085e+01, 5.927e+00,
#     5.000e+00, 3.950e+00, 3.850e+00, 3.825e+00, 3.850e+00,
#     3.900e+00, 3.975e+00, 4.065e+00, 4.200e+00, 4.300e+00,
#     4.425e+00, 4.575e+00, 4.725e+00, 4.825e+00, 4.900e+00,
#     4.950e+00, 5.025e+00, 5.150e+00, 5.225e+00, 5.250e+00,
#     5.225e+00, 5.100e+00, 4.750e+00, 4.200e+00, 3.500e+00,
#     2.825e+00, 2.050e+00, 1.330e+00, 8.500e-01, 5.400e-01,
#     4.000e-01, 3.400e-01, 2.800e-01, 2.400e-01, 2.000e-01,
#     ])*1.0e-6

# Mid-latitude summer
p_std = array([
    1.013E+03, 9.020E+02, 8.020E+02, 7.100E+02, 6.280E+02,
    5.540E+02, 4.870E+02, 4.260E+02, 3.720E+02, 3.240E+02,
    2.810E+02, 2.430E+02, 2.090E+02, 1.790E+02, 1.530E+02,
    1.300E+02, 1.110E+02, 9.500E+01, 8.120E+01, 6.950E+01,
    5.950E+01, 5.100E+01, 4.370E+01, 3.760E+01, 3.220E+01,
    2.770E+01, 1.907E+01, 1.320E+01, 9.300E+00, 6.520E+00,
    4.640E+00, 3.330E+00, 2.410E+00, 1.760E+00, 1.290E+00,
    9.510E-01, 5.150E-01, 2.720E-01, 1.390E-01, 6.700E-02,
    3.000E-02, 1.200E-02, 4.480E-03, 1.640E-03, 6.250E-04,
    2.580E-04, 1.170E-04, 6.110E-05, 3.560E-05, 2.270E-05,
    ])*1e+2
h2o_vmr_std = array([
    1.876E+04, 1.378E+04, 9.680E+03, 5.984E+03, 3.813E+03,
    2.225E+03, 1.510E+03, 1.020E+03, 6.464E+02, 4.129E+02,
    2.472E+02, 9.556E+01, 2.944E+01, 8.000E+00, 5.000E+00,
    3.400E+00, 3.300E+00, 3.200E+00, 3.150E+00, 3.200E+00,
    3.300E+00, 3.450E+00, 3.600E+00, 3.850E+00, 4.000E+00,
    4.200E+00, 4.450E+00, 4.700E+00, 4.850E+00, 4.950E+00,
    5.000E+00, 5.100E+00, 5.300E+00, 5.450E+00, 5.500E+00,
    5.500E+00, 5.350E+00, 5.000E+00, 4.400E+00, 3.700E+00,
    2.950E+00, 2.100E+00, 1.330E+00, 8.500E-01, 5.400E-01,
    4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01,
    ])*1e-6

################################################################################
# Function to calculate mass mixing ratios
################################################################################

def check_flag(flag):
  if flag not in mean_molar_weight.keys():
    print('Error: Flag ' + flag + ' not supported.')
    exit()

# Calculates mass mixing ratio
def mass_mixing_ratio(gas, P, T, flag='', return_dry = False,
    mix_ratio_override = {}):

  # Check if flag is supported
  check_flag(flag)

  # Set number mixing ratio

  if gas == 'H2O':
    ip = interp1d(log10(p_std), log10(h2o_vmr_std), 
        bounds_error = False, fill_value = 'extrapolate')
    mmr = 10.0**ip(log10(P))

  elif gas == 'CO2':
    if flag == 'setup_1' or flag == 'setup_2' or flag == 'setup_3':
      mmr = 0.01*ones(len(P))
    elif flag == 'setup_4':
      mmr = 0.89*ones(len(P))

  elif gas == 'CH4':
    mmr = 0.01*ones(len(P))

  elif gas == 'N2':
    if flag == 'setup_1':
      mmr = 0.98*ones(len(P))
    elif flag == 'setup_3':
      mmr = 0.49*ones(len(P))
    else:
      mmr = 0.0*ones(len(P))

  elif gas == 'H2':
    if flag == 'setup_2':
      mmr = 0.98*ones(len(P))
    elif flag == 'setup_3':
      mmr = 0.49*ones(len(P))
    elif flag == 'setup_4':
      mmr = mmr = 0.1*ones(len(P))
    else:
      mmr = 0.0*ones(len(P))

  else:
    return 0.0*ones(len(P))

  # Calculate mass mixing ratio from number mixing ratio
  mmr = mmr*molar_weight[gas]/mean_molar_weight[flag]

  # Convert the dry mass mixing ratio to specific if requested
  if not return_dry:
    mmr = mmr/(1 + mass_mixing_ratio('H2O', P, T, flag = flag,
        return_dry = True, mix_ratio_override = mix_ratio_override))

  return mmr

# Returns maximum allowed mixing ratios for each species
def max_mass_mixing_ratio(gas, flag=flag, return_dry = False,
    mix_ratio_override = {}):

  max_mass_mixing_ratio = mass_mixing_ratio(gas, 0., 0., flag = flag,
      return_dry = return_dry, mix_ratio_override = mix_ratio_override)

  return max_mass_mixing_ratio
