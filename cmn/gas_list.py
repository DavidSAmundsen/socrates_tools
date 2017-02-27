from pylab import *

# Mean molecular weights in grams per mole from NIST
molar_weight = {
    'H':        1.00794,
    'H2':       2.01588,
    'He':       4.002602,
    'H2O':      18.0153,
    'CO':       28.0106,
    'CO2':      44.0095,
    'CH4':      16.0430,
    'O2':       31.9988,
    'NH3':      17.0306,
    'TiO':      63.866,
    'VO':       66.9409,
    'Na':       22.98976928,
    'K':        39.0983,
    'FeH':      56.853,
    'CrH':      53.004,
    'Li':        6.941,
    'Rb':       85.4678,
    'Cs':      132.9054519,
    'PH3':      33.99758,
    'C2H2':     26.0373,
    'HCN':      27.0253,
    'O3':       47.9982,
    'N2O':      44.0128,
    'H2S':      34.081,
    'AIR':      28.966,
    }

# Long names for gases
gas_name = {
    'H2':   'Hydrogen',
    'He':   'Helium',
    'H2O':  'Specific humidity',
    'CO':   'Carbon monoxide',
    'CO2':  'Carbon dioxide',
    'CH4':  'Methane',
    'O2':   'Oxygen',
    'NH3':  'Ammonia',
    'TiO':  'Titanium oxide',
    'VO':   'Vanadium oxide',
    'Na':   'Sodium',
    'K':    'Potassium',
    'FeH':  'Iron hydride',
    'CrH':  'Chromium hydride',
    'Li':   'Lithium',
    'Rb':   'Rubidium',
    'Cs':   'Cesium',
    'PH3':  'Phosphine',
    'C2H2': 'Acetylene',
    'HCN':  'Hydrogen cyanide',
    'O3':   'Ozone',
    'N2O':  'Nitrous oxide',
    'H2S':  'Hydrogen sulphide',
    'AIR':  'Dry air',
    }

# HITRAN molecule IDs
gas_name_hitran = array([
    None,
    'H2O',
    'CO2',
    'O3',
    'N2O',
    'CO',
    'CH4',
    'O2',
    ])

# Get HITRAN ID for gases
def gas_id_hitran(gas_name):

# Convert to array if only a string
  if (type(gas_name) == str):
    gas_name = array([gas_name])

# Loop through all gas names and find HITRAN ID
  gas_id = zeros(len(gas_name), dtype=int)
  for i in arange(len(gas_name)):
    gas_id[i] = argwhere(gas_name_hitran == gas_name[i])[0,0]
  
# Only return an integer if only a gas was provided
  if (len(gas_id) == 1):
    return gas_id[0]
  else:
    return gas_id
