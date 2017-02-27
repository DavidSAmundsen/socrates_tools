# Contains a function  to calculate mass mixing ratios for solar abundance
# atmospheres such as hot Jupiters and brown dwarfs. Uses the abundances from
# Burrows & Sharp 1999 with parameterisations added for additional gases.

from pylab import *

# External files with data required
from constants import r_gas_cal
from gas_list import molar_weight

# Elemental number fractions of hydrogen and helium
elem_num_frac_h = 0.91183
elem_num_frac_he = 1.0 - elem_num_frac_h

# Elemental abundances in units of log10(N_X/N/H) + 12
elem_abund = {
    'H':  12.00,
    'He': 10.98,
    'Li':  3.26,
    'C':   8.50,
    'N':   7.86,
    'O':   8.76,
    'Na':  6.24,
    'Si':  7.51,
    'P':   5.46,
    'S':   7.16,
    'K':   5.03,
    'Ti':  4.95,
    'V':   3.93,
    'Cr':  5.65,
    'Fe':  7.52,
    'Rb':  2.52,
    'Cs':  1.08,
    }

# Calculate H2 and He number fractions
num_frac_h2 = (elem_num_frac_h/2.0)/(elem_num_frac_h/2.0 + elem_num_frac_he)
num_frac_he = 1.0 - num_frac_h2

# Mean molecular weight of atmosphere
mean_molar_weight = (molar_weight['H2']*num_frac_h2 +
    molar_weight['He']*num_frac_he)

# Calculate number fractions from elemental abundances relative to H
elem_num_frac = dict(elem_abund)
for key in elem_abund.keys():
  elem_num_frac[key] = 10**(elem_abund[key] - 12.0)

# Chemical equilibrium constants
eqcoeff_1 = array([
     1.106131e+6,
    -5.6895e+4,
    62.565e+0,
    -5.81396e-4,
     2.346515e-8
    ])
eqcoeff_2 = array([
     8.16413e+5,
    -2.9109e+4,
    58.5878e+0,
    -7.8284e-4,
     4.729048e-8
    ])

# Temperature below which oxygen is removed by silicon atoms
T_remove_O = 1500.0

# Average number of oxygen atoms removed per silicon atom
x_Si = 3.28

# Transition scales for condensing species
T_trans_scale = {
    'O':   10.0,
    'TiO': 10.0,
    'VO':  10.0,
    'Na':  20.0,
    'K':   20.0,
    'Li':  20.0,
    'Rb':  20.0,
    'Cs':  20.0,
    }

# Condensation/transofmration curve parameterisations
poly_highp = {
    'TiO': array([-3.96274342e-05, 5.20741797e-04]),
    'VO':  array([-3.96274342e-05, 5.20741797e-04]),
    'Na':  array([-5.58256919e-05, 8.81851644e-04]),
    'K':   array([-5.46977180e-05, 8.19104478e-04]),
    'Li':  array([-3.50995394e-05, 6.51993843e-04]),
    'Rb':  array([-6.06654087e-05, 8.09569974e-04]),
    'Cs':  array([-5.29210264e-05, 7.71577097e-04]),
}
poly_lowp = {
    'TiO': array([-2.91788901e-05, 5.11801742e-04]),
    'VO':  array([-2.91788901e-05, 5.11801742e-04]),
    'Na':  array([-6.69921629e-05, 8.90116829e-04]),
    'K':   array([-6.46633572e-05, 8.29549449e-04]),
    'Li':  array([-3.55469185e-05, 6.52116945e-04]),
    'Rb':  array([-3.19328287e-05, 8.69542964e-04]),
    'Cs':  array([-3.85306167e-05, 7.63040762e-04]),
    }
p_bar_trans = {
    'TiO': 10.0,
    'VO':  10.0,
    'Na':  10.0,
    'K':   10.0,
    'Li':  10.0,
    'Rb':  1.0e-2,
    'Cs':  10.0,
    }

################################################################################
# Function to calculate mass mixing ratios
################################################################################

# Returns the temperature given a polynomial fit
def cond_curve(log_p_bar, poly):
  return 1.0/polyval(poly, log_p_bar)

# Calculates condensation temperature of TiO and VO
def cond_temp(P, gas):
  log_p_bar = log10(P/1.0e+5)
  work = 1.0/(exp(-(log_p_bar - log10(p_bar_trans[gas]))/0.1) + 1.0)
  t_cond = cond_curve(log_p_bar, poly_lowp[gas])*(1.0 - work) + \
      cond_curve(log_p_bar, poly_highp[gas])*work

  return t_cond

# Calculates mass mixing ratio
def mass_mixing_ratio(gas, P, T, flag='', return_dry=True):

# Calculate gas partial pressure relative to H2 partial pressure

  if (gas == 'H2'):
      mass_mixing_ratio = num_frac_h2/num_frac_h2

  elif (gas == 'He'):
      mass_mixing_ratio = num_frac_he/num_frac_h2

  elif (gas == 'CO' or gas == 'CH4' or gas == 'H2O'):

    eq_cnst = exp((eqcoeff_1[0]/T + eqcoeff_1[1] + eqcoeff_1[2]*T +
        eqcoeff_1[3]*T**2 + eqcoeff_1[4]*T**3)/(r_gas_cal*T))

    frac_remove_O = (1.0e+0/
        (exp((T - T_remove_O)/T_trans_scale['O']) + 1.0e+0))

    elem_num_frac_O_depleted = (elem_num_frac['O'] -
        x_Si*elem_num_frac['Si']*frac_remove_O)

    work = (elem_num_frac['C'] + elem_num_frac_O_depleted +
        (num_frac_h2*P/101325.0)**2/(2.0*eq_cnst))

    mass_mixing_ratio = work - sqrt(work**2.0 -
        4.0*elem_num_frac['C']*elem_num_frac_O_depleted)

    if (gas == 'CH4'):
      mass_mixing_ratio = 2.0*elem_num_frac['C'] - mass_mixing_ratio

    elif (gas == 'H2O'):
      mass_mixing_ratio = 2.0*elem_num_frac_O_depleted - mass_mixing_ratio

  elif (gas == 'NH3'):

    eq_cnst = exp((eqcoeff_2[0]/T + eqcoeff_2[1] + eqcoeff_2[2]*T +
        eqcoeff_2[3]*T**2 + eqcoeff_2[4]*T**3)/(r_gas_cal*T))

    work = (num_frac_h2*P/101325.0)**2.0/(8.0*eq_cnst)

    mass_mixing_ratio = 2.0*(sqrt(work*(2.0*elem_num_frac['N'] +
        work)) - work)

  elif (gas == 'TiO' ):
    mass_mixing_ratio = 2.0*elem_num_frac['Ti']/(
        exp(-(T - cond_temp(P, gas))/T_trans_scale[gas]) + 1.0)

  elif (gas == 'VO'):
    mass_mixing_ratio = 2.0*elem_num_frac['V']/(
        exp(-(T - cond_temp(P, gas))/T_trans_scale[gas]) + 1.0)

  elif (gas == 'Na' or gas == 'K' or gas == 'Li' or gas == 'Rb' or gas == 'Cs'):
    mass_mixing_ratio = 2.0*elem_num_frac[gas]/(
        exp(-(T - cond_temp(P, gas))/T_trans_scale[gas]) + 1.0)

  else:
    print('Warning: Gas ' + gas + ' not found in mass_mixing_ratio.')
    return 0.0

# Calculate mass mixing ratio
  mass_mixing_ratio = (mass_mixing_ratio*num_frac_h2*
      molar_weight[gas]/mean_molar_weight)

  return mass_mixing_ratio


# Returns maximum allowed mixing ratios for each species
def max_mass_mixing_ratio(gas):

# Set maximum gas partial pressure relative to H2 partial pressure
  if (gas == 'H2'):
    max_mass_mixing_ratio = 1.0
  elif (gas == 'He'):
    max_mass_mixing_ratio = num_frac_he/num_frac_h2
  elif (gas == 'CO' or gas == 'CH4' or gas == 'CO2'):
    max_mass_mixing_ratio = 2.0*elem_num_frac['C']
  elif (gas == 'H2O'):
    max_mass_mixing_ratio = 2.0*elem_num_frac['O']
  elif (gas == 'NH3'):
    max_mass_mixing_ratio = 2.0*elem_num_frac['N']
  elif (gas == 'TiO'):
    max_mass_mixing_ratio = 2.0*elem_num_frac['Ti']
  elif (gas == 'VO'):
    max_mass_mixing_ratio = 2.0*elem_num_frac['V']
  elif (gas == 'Na' or gas == 'K' or gas == 'Li' or
      gas == 'Rb' or gas == 'Cs'):
    max_mass_mixing_ratio = 2.0*elem_num_frac[gas]
  elif (gas == 'FeH'):
    max_mass_mixing_ratio = 2.0*elem_num_frac['Fe']
  elif (gas == 'CrH'):
    max_mass_mixing_ratio = 2.0*elem_num_frac['Cr']
  elif (gas == 'PH3'):
    max_mass_mixing_ratio = 2.0*elem_num_frac['P']
  elif (gas == 'C2H2'):
    max_mass_mixing_ratio = 7.0e-3/num_frac_h2 # From Pascal
  elif (gas == 'HCN'):
    max_mass_mixing_ratio = 3.0e-3/num_frac_h2 # From Pascal
  elif (gas == 'H2S'):
    max_mass_mixing_ratio = 2.0*elem_num_frac['S']
  else:
    print('Error: Gas ' + gas + ' not found in max_mass_mixing_ratio.')
    return 0.0

# Calculate mass mixing ratio
  max_mass_mixing_ratio = (max_mass_mixing_ratio*num_frac_h2*
      molar_weight[gas]/mean_molar_weight)

  return max_mass_mixing_ratio
