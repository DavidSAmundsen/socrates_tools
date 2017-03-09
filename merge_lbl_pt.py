# This script merges line-by-line files for different P-T grids.

from pylab import *
from netCDF4 import *
import os
import sys

lbl_dir = '/scr2/socrates/data/abs_coeff'
file_in = [
    'co2_sblbl_lw_pt663_gf1.0.nc',
    'co2_sblbl_lw_pt130_gf1.0.nc',
    ]
file_out = 'co2_sblbl_lw_pt793_gf1.0.nc'
# file_in = [
#     'h2o_lbl_lwf_pt663.nc',
#     'h2o_lbl_lwf_pt130.nc',
#     ]
# file_out = 'h2o_lbl_lwf_pt793.nc'

# Get combined P-T grid
p_calc = array([])
t_calc = array([])
for i in arange(len(file_in)):
   fid_in = Dataset(os.path.join(lbl_dir, file_in[i]), 'r')
   p_calc = append(p_calc, fid_in.variables['p_calc'][:])
   t_calc = append(t_calc, fid_in.variables['t_calc'][:])
   fid_in.close()

# Use the first file to get wavenumber and potentially gas fractions points
fid_in = Dataset(os.path.join(lbl_dir, file_in[0]), 'r')
nu = fid_in.variables['nu'][:]
nu_step = fid_in.variables['nu'].step
if 'gas_frac' in fid_in.variables.keys():
  gas_frac = fid_in.variables['gas_frac'][:]
fid_in.close()

# Add everything except opacities to output file
fid_out = Dataset(os.path.join(lbl_dir, file_out), 'w')
fid_out.createDimension('nu', len(nu))
fid_out.createDimension('pt_pair', len(p_calc))
p_calc_out = fid_out.createVariable('p_calc', 'f4', ('pt_pair'))
t_calc_out = fid_out.createVariable('t_calc', 'f4', ('pt_pair'))
nu_out = fid_out.createVariable('nu', 'f8', ('nu'))
p_calc_out[:] = p_calc
t_calc_out[:] = t_calc
nu_out[:] = nu
if 'gas_frac' in locals():
  fid_out.createDimension('gas_frac', len(gas_frac))
  gas_frac_out = fid_out.createVariable('gas_frac', 'f4', ('gas_frac'))
  kabs_out = fid_out.createVariable('kabs', 'f4', ('gas_frac', 'pt_pair', 'nu'))
  gas_frac_out[:] = gas_frac
else:
  kabs_out = fid_out.createVariable('kabs', 'f4', ('pt_pair', 'nu'))

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
kabs_out.title = 'absorption'
kabs_out.long_name = 'absorption'
kabs_out.units = 'm2 kg-1'
if 'gas_frac' in locals():
  gas_frac_out.title = 'gas fraction'
  gas_frac_out.long_name = 'gas fraction'

# Loop through input files and add opacities for each P-T point
i_pt = -1
for i in arange(len(file_in)):
  print('Processing ' + file_in[i] + '...')
  fid_in = Dataset(os.path.join(lbl_dir, file_in[i]), 'r')

  # Check wavenumber and gas fractions
  if any(fid_in.variables['nu'][:] != nu):
    print('Wavenumbers do not match.')
    sys.exit(1)
  elif 'gas_frac' in locals():
    if any(fid_in.variables['gas_frac'][:] != gas_frac):
      print('Gas fractions do not match.')
      sys.exit(1)

  for j in arange(fid_in.dimensions['pt_pair'].size):
    i_pt = i_pt + 1
    if 'gas_frac' in locals():
      for k in arange(len(gas_frac)):
        kabs_out[k,i_pt,:] = fid_in.variables['kabs'][k,j,:]
    else:
      kabs_out[i_pt,:] = fid_in.variables['kabs'][j,:]
  fid_in.close()
fid_out.close()
print('Done.')