# This script merges line-by-line files for different gas fractions.

from pylab import *
from netCDF4 import *
import os
import sys

lbl_dir = '/scr2/socrates/data/abs_coeff'
file_in = [
    'co2_sblbl_lw_pt663_gf0.1.nc',
    'co2_sblbl_lw_pt663_gf1.0.nc',
    ]
file_out = 'co2_sblbl_lw_pt663.nc'

# Get combined gas fractions
gas_frac = array([])
for i in arange(len(file_in)):
   fid_in = Dataset(os.path.join(lbl_dir, file_in[i]), 'r')
   gas_frac = append(gas_frac, fid_in.variables['gas_frac'][:])
   fid_in.close()

# Use the first file to get wavenumber and P-T points
fid_in = Dataset(os.path.join(lbl_dir, file_in[0]), 'r')
p_calc = fid_in.variables['p_calc'][:]
t_calc = fid_in.variables['t_calc'][:]
nu = fid_in.variables['nu'][:]
nu_step = fid_in.variables['nu'].step
fid_in.close()

# Add everything except opacities to output file
fid_out = Dataset(os.path.join(lbl_dir, file_out), 'w')
fid_out.createDimension('nu', len(nu))
fid_out.createDimension('pt_pair', len(p_calc))
fid_out.createDimension('gas_frac', len(gas_frac))
p_calc_out = fid_out.createVariable('p_calc', 'f4', ('pt_pair'))
t_calc_out = fid_out.createVariable('t_calc', 'f4', ('pt_pair'))
nu_out = fid_out.createVariable('nu', 'f8', ('nu'))
gas_frac_out = fid_out.createVariable('gas_frac', 'f4', ('gas_frac'))
kabs_out = fid_out.createVariable('kabs', 'f4', ('gas_frac', 'pt_pair', 'nu'))
p_calc_out[:] = p_calc
t_calc_out[:] = t_calc
nu_out[:] = nu
gas_frac_out[:] = gas_frac

# Add attributes
p_calc_out.title = 'pressure'
p_calc_out.long_name = 'pressure'
p_calc_out.units = 'Pa'
t_calc_out.title = 'temperature'
t_calc_out.long_name = 'temperature'
t_calc_out.units = 'K'
nu_out.title = 'wavenumber'
nu_out.long_name = 'wavenumber'
nu_out.units = 'm-1'
nu_out.step = nu_step
gas_frac_out.title = 'gas fraction'
gas_frac_out.long_name = 'gas fraction'
kabs_out.title = 'absorption'
kabs_out.long_name = 'absorption'
kabs_out.units = 'm2 kg-1'

# Loop through input files and add opacities for each gas fraction
i_gf = -1
for i in arange(len(file_in)):
  print('Processing ' + file_in[i])
  fid_in = Dataset(os.path.join(lbl_dir, file_in[i]), 'r')

  # Check wavenumber and P-T points
  if any(fid_in.variables['nu'][:] != nu):
    print('Wavenumbers do not match.')
    sys.exit(1)
  elif any(fid_in.variables['p_calc'][:] != p_calc):
    print('Pressures do not match.')
    sys.exit(1)
  elif any(fid_in.variables['t_calc'][:] != t_calc):
    print('Temperatures do not match.')
    sys.exit(1)

  for j in arange(fid_in.dimensions['gas_frac'].size):
    i_gf = i_gf + 1
    for k in arange(len(p_calc)):
      kabs_out[i_gf,k,:] = fid_in.variables['kabs'][j,k,:]
  fid_in.close()
fid_out.close()
print('Done.')