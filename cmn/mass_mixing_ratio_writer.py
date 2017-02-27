from pylab import *
from scipy.interpolate import interp1d
from sys import exit
from constants import n_loschmidt, n_avogadro

################################################################################
# Data
################################################################################

# External files with data required
from gas_list import molar_weight

mean_molar_weight = 28.964
grav_acc = 9.80665

# Column indices for ULGAS data
i_gas = {
    'H2O':    1,
    'CO2':    2,
    'O3':     3,
    'N2O':    4,
    'CH4':    5,
    'CFC-11': 6,
    'CFC-12': 7,
    'NO2':    8,
    }

################################################################################
# Function to calculate mass mixing ratios
################################################################################

# Calculates mass mixing ratio
def mass_mixing_ratio(gas, P, T, flag='', return_dry = False,
    mix_ratio_override = {}):
  
  # Set dry mass mixing ratios
  if (gas in mix_ratio_override.keys()):
    # Override mass mixing ratio
    mmr = mix_ratio_override[gas]
  else:
    mmr = (vol_mix_ratio(gas, P, T, flag = flag)*
        molar_weight[gas]/mean_molar_weight)

  # Convert to specific mixing ratio if requested
  if (not return_dry):
    mmr = mmr/(1 + mass_mixing_ratio('H2O', P, T, flag = flag,
        return_dry = True, mix_ratio_override = mix_ratio_override))

  return mmr
  
def vol_mix_ratio(gas, P, T, flag = ''):

  # No file with abundances if no flag
  if not flag:
    sys.exit('No WRITER file with mixing ratios provided')

  # Read file with ULGAS data
  data = loadtxt(flag, skiprows=2, usecols=(1,7,8,9,10,11,12,13,14))

  # Read file with pressure at levels data
  plev = loadtxt(flag[:-10] + flag[-4:], skiprows=2, usecols=(1,))

  # Convert pressure to Pa
  data[:,0] = data[:,0]*1.0e+2
  plev = plev*1.0e+2

  # Convert cm at STP to volume mixing ratios
  for i in arange(1,size(data,1)):
    data[:,i] = (data[:,i]*1.0e-5*n_loschmidt*mean_molar_weight/
        (n_avogadro*diff(plev)/grav_acc))

  # If gas is in WRITER table, interpolate to correct pressure
  if (gas in i_gas.keys()):
    vmr_func = interp1d(log10(data[:,0]), data[:,i_gas[gas]],
        bounds_error = False, fill_value = 'extrapolate')
    vmr = vmr_func(log10(P))
    
    # Make sure all mixing ratios are positive
    vmr[argwhere(vmr < 0.0)] = 0.0
    
    return vmr
    
  # Otherwise set mixing ratio to zero
  else:
    return 0.0

  return 0.0