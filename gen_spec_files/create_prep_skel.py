from pylab import *
from netCDF4 import *
import os
import pickle
import sys

################################################################################
# Modern Earth, for short-wave spectral file
################################################################################

absorbers = array([1, 2, 3, 4, 6, 7, 9, 25])
continua = array([1])
aerosols = array([])
sp_id = 'dsa' # Identifier for spectral file
sp_region = 'sw'

grid_type = 'const_dlambda_lambda' # Grid of bands to be used
bands_min = 2.0e-7 # Minimum wavelength in m
bands_max = 2.0e-5 # Maximum wavelength m
n_bands = 500 # Number of bands

gases_in_bands = False # Do not include a list of gases in each band

################################################################################
# Archean Earth with H2O, CO2 and CH4
################################################################################

# nu_precision = resolution/10.0 # Precision of gases in bands, should be < resolution
# absorbers = array([1, 2, 6])
# continua = array([1])
# aerosols = array([])
# sp_id = 'ar_dsa' # Identifier for spectral file
# sp_region = 'lw'
# bands_min = 0.0 # resolution/10.0 # Minimum wavenumber in cm-1
# bands_max = 3500.0 # Maximum wavenumber cm-1
# gases_in_bands = False # Do not include a list of gases in each band
# 
# path_abs_coeff = '/home/damundse/Abs_coeff/'
# absorber_names = ['H2O', 'CO2', 'CH4']
# note_abs = ''

################################################################################

# Maximum wavenumber and minimum wavelength for which to include water continuum
max_wn_h2o_continuum = 20000.
min_wl_h2o_continuum = 1.0e-2/max_wn_h2o_continuum

# Construct band limits
if grid_type == 'const_nu':
  band_min = arange(bands_min, bands_max, resolution)
  band_max = arange(bands_min+resolution, bands_max+resolution, resolution)

  # Have band 1 start at 1.0 cm-1 as 0.0 gives infinite wavelength
  band_min[0] = 1.0

  # Calculate number of bands
  n_bands = len(band_min)

  # Find the bands in which to include water continuum
  l_include_continuum = band_max <= max_wn_h2o_continuum

elif grid_type == 'const_dlambda_lambda':
  bands = logspace(log10(bands_min), log10(bands_max), n_bands + 1)
  band_min = bands[:-1]
  band_max = bands[1:]

  # Band limits should be on integer wavenumbers in m:
  band_min = 1.0/around(1.0/band_min)
  band_max = 1.0/around(1.0/band_max)

  # If spectral region is short-wave, then have a clear division between visible
  # and near-IR at min_wl_h2o_continuum
  if sp_region == 'sw':
    # Find band limit closest to max_wn_h2o_continuum
    i_vis_nir = argmin(abs(band_min - min_wl_h2o_continuum))
    
    # Adjust band limit i_vis_nir to match min_wl_h2o_continuum
    band_min[i_vis_nir] = min_wl_h2o_continuum
    band_max[i_vis_nir-1] = min_wl_h2o_continuum

  # Find the bands in which to include water continuum
  l_include_continuum = band_min >= min_wl_h2o_continuum
else:
  print('Unrecognised grid type.')
  sys.exit(1)

# Construct spectral file name
sp_name = 'sp_' + sp_region + '_{:g}_'.format(n_bands) + sp_id + '_skel'
mk_sp_name = 'mk_' + sp_name

# Read file with previously generated data
if gases_in_bands:
  sv_file = 'prep_skel_data/' + mk_sp_name + '.pickle'
  if os.path.isfile(sv_file):
    with open('prep_skel_data/' + mk_sp_name + '.pickle', 'rb') as f:
      (absorber_names_sv, sources_sv, bands_min_sv, bands_max_sv, resolution_sv,
          absorbers_nu_min_sv, absorbers_nu_max_sv) = pickle.load(f)
        
    # Use the data in the file only if the band limits used match
    if (bands_min == bands_min_sv and bands_max == bands_max_sv and
        resolution == resolution_sv):
      use_sv = True
    else:
      use_sv = False

  else:
    use_sv = False

  # Find for each absorber the min and max wavenumber where data exists
  absorbers_nu_min = zeros(len(absorbers))
  absorbers_nu_max = zeros(len(absorbers))
  print('Absorber, nu_min, nu_max')
  for i in range(len(absorbers)):
    if use_sv:
      ind_sv = argwhere(absorber_names[i] == array(absorber_names_sv))
      if ind_sv.size:
        use_sv_i = True
      else:
        use_sv_i = False
    else:
      use_sv_i = False
  
    if use_sv_i:
      ind_sv = ind_sv[0,0]
      absorbers_nu_min[i] = absorbers_nu_min_sv[ind_sv]
      absorbers_nu_max[i] = absorbers_nu_max_sv[ind_sv]
    else:
      if ('sources' in locals()):
        fid = Dataset(path_abs_coeff + 'abs_coeff_' + absorber_names[i] + '_' +
            sources[i] + note_abs + '.nc', 'r')
      else:
        fid = Dataset(path_abs_coeff + absorber_names[i].lower() +
            '_lbl_lw.nc', 'r')
      i_step = int(round(nu_precision/(fid.variables['nu'].step*1e-2)))
      nu = fid.variables['nu'][::i_step]*1e-2
      kabs = fid.variables['kabs'][-1,::i_step]

      absorbers_nu_min[i] = max(nu[min(argwhere(kabs > 0.))] - nu_precision,
          band_min[0])
      absorbers_nu_max[i] = min(nu[max(argwhere(kabs > 0.))] + nu_precision,
          band_max[-1])
      fid.close()
  
    print(absorber_names[i], absorbers_nu_min[i], absorbers_nu_max[i])

print('Constructing ' + mk_sp_name + ' file...')

# Open prep_spec script and write header
fid = open(mk_sp_name, 'w')
fid.write('# This script creates a basic skeleton {:g} band '.format(n_bands) +
    'spectral file.\n\n')
fid.write('rm -f ' + sp_name + '\n\n')

# Start writing call to prep_spec
fid.write('prep_spec > /dev/null <<EOF\n')
fid.write(sp_name + '\n')
fid.write('{:g}\n'.format(n_bands))

# Number of absorbers, list of absorber indices and number of aerosols
fid.write('{:g}\n'.format(absorbers.size))
if (absorbers.size):
  for i in range(absorbers.size):
    fid.write('{:g}\n'.format(absorbers[i]))
fid.write('{:g}\n'.format(aerosols.size))
if (aerosols.size):
  for i in range(aerosols.size):
    fid.write('{:g}\n'.format(aerosols[i]))

# Set unit of band limits
if grid_type == 'const_nu':
  fid.write('c\n')
elif grid_type == 'const_dlambda_lambda':
  fid.write('m\n')

# Band limits
for i in range(len(band_min)):
  fid.write('{:1.9E}     {:1.9E}\n'.format(band_max[i],
      band_min[i]).replace('E', 'D'))

# Gases present in each band
absorber_band_min = ones(len(absorbers), dtype=int)*-1
absorber_band_max = ones(len(absorbers), dtype=int)*-1
for i in range(n_bands):
  if gases_in_bands:
    n_abs_band = 0
    for j in range(len(absorbers)):
      if ((absorbers_nu_min[j] <= band_min[i] and
          absorbers_nu_max[j] >= band_max[i])):
        if (n_abs_band > 0):
          fid.write(' ')
        fid.write('{:g}'.format(absorbers[j]))
        n_abs_band = n_abs_band + 1
      
        absorber_band_max[j] = i
        if (absorber_band_min[j] < 0):
          absorber_band_min[j] = i
  else:
    fid.write('{:g}'.format(0))
  fid.write('\n')

# Change band indexing to start on 1
absorber_band_min = absorber_band_min + 1
absorber_band_max = absorber_band_max + 1
  
# Type of continua in bands
for i in range(n_bands):
  if (size(continua) > 0 and l_include_continuum[i]):
    # There are continuum absorbers
    n_cont_band = 0
    for j in range(len(continua)):
      if (n_cont_band > 0):
        fid.write(' ')
      fid.write('{:g}'.format(continua[j]))
      n_cont_band = n_cont_band + 1
    fid.write('\n')
  else: # No continuum absorbers
    fid.write('0\n')

# Do not exclude regions from bands
fid.write('n\n')

# Finish
fid.write('-1\n')
fid.write('EOF')

fid.close()

print('Done.')

if gases_in_bands:
  print('Absorber, First band, Last band')
  for i in range(len(absorbers)):
    print(absorber_names[i], absorber_band_min[i], absorber_band_max[i])

  # Crate updated pickle files if number of gases has been increased
  if (not use_sv or
      len(absorber_names) > len(absorber_names_sv)):
    with open('prep_skel_data/' + mk_sp_name + '.pickle', 'wb') as f:
      if ('sources' not in locals()):
        sources = ['HITRAN']*len(absorber_names)
      pickle.dump([
          absorber_names, sources, bands_min, bands_max, resolution,
          absorbers_nu_min, absorbers_nu_max,
          ], f)
