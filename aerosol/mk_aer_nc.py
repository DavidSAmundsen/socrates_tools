from pylab import *
from netCDF4 import *
import os
import glob
from aer_density_component import aer_density_component

################################################################################
# User definitions
################################################################################

spectral_file_name = 'sp_sw_21_dsa'
spectral_file_folder = 'sp_sw_dsa'
solar_spec='sun'

################################################################################
# Diagnostics, paths and aerosol types
################################################################################

# Spectral file for diagnostics
spectral_file_name_550nm = 'sp_diag'

socrates_dir = '/scr2/socrates'
scatter_base_dir = os.path.join(socrates_dir, 'data', 'aerosol')
spectral_file_dir = os.path.join(socrates_dir, 'spectral_files')
if spectral_file_name[3:5] == 'lw':
  scatter_data_dir = os.path.join(scatter_base_dir, spectral_file_name)
  aer_nc = os.path.join(spectral_file_dir, spectral_file_folder,
      'aer_' + spectral_file_name[3:] + '.nc')
elif spectral_file_name[3:5] == 'sw':
  scatter_data_dir = os.path.join(scatter_base_dir, spectral_file_name,
      solar_spec)
  aer_nc = os.path.join(spectral_file_dir, spectral_file_folder,
      'aer_' + spectral_file_name[3:] + '_' + solar_spec + '.nc')
else:
  print('Spectral region not recognised.')
  sys.exit(1)
scatter_data_550nm_dir = os.path.join(scatter_base_dir,
    spectral_file_name_550nm)

# List of aerosol types
aer_type = array([
    'sdust',
    'sinyukdust',
    'h2so4',
    'soot',
    'sulphate',
    'seasalt',
    'antwater',
    'organic',
    ])

################################################################################
# Functions
################################################################################

def get_radius(scatter_data_dir, aer_type):

  # Get list of all files with data for aer_type
  aer_file = glob.glob(os.path.join(scatter_data_dir, aer_type) + '_*')

  # Create list of radii
  radius_str = [0]*len(aer_file)
  radius = zeros(len(aer_file))
  for i_file in arange(len(aer_file)):
    radius_str[i_file] = aer_file[i_file].split('/')[-1].split('_')[-1][:-4]
    radius[i_file] = float(radius_str[i_file])
  radius_str = array(radius_str)

  # Sort radii in ascending order
  i_sorted = argsort(radius)
  radius_str = radius_str[i_sorted]

  return radius_str


def get_n_band(scatter_file, l_humidity):

  fid = open(scatter_file)

  # Skip to data
  while True:
    line = fid.readline()
    if line[:10] == '*FILE TYPE':
      break

  # Skip until start of band-dependent data
  if l_humidity:
    for i in arange(17): fid.readline()
  else:
    for i in arange(13): fid.readline()

  # Every other line is a new band until end of data
  n_band = 0
  n_newline = 0
  while True:
    line = fid.readline()
    if line[:10] == '*FILE TYPE' or n_newline > 1:
      break
    elif line == '\n':
      n_newline = n_newline + 1
    else:
      n_band = n_band + 1
      n_newline = 0

  return n_band


def humidity_dependence(fname):

  fid = open(fname, 'r')

  # Look for file type
  while True:
    line = fid.readline()
    if line[:18] == '*FILE TYPE =     8':
      # No humidity dependence
      fid.close()
      return False, 1
    elif line[:18] == '*FILE TYPE =     9':
      # Data has humidity dependnece
      break

  # Count number of humidities in file
  n_humidity = 1
  while True:
    line = fid.readline()
    if line[:18] == '*FILE TYPE =     9':
      # Another humidity
      n_humidity = n_humidity + 1
    elif line == '':
      # End of file
      break

  fid.close()
  return True, n_humidity


def read_scatter_data(fid, l_humidity):

  # Skip to data
  while True:
    line = fid.readline()
    if line[:10] == '*FILE TYPE':
      break

  # Skip data not needed
  for i in arange(5): fid.readline()

  # Get aerosol component and name
  line = fid.readline()
  aer_cmp = int(line[24:29])
  aer_name = line[33:53]

  # Read humidity
  if l_humidity:
    fid.readline()
    line = fid.readline()
    humidity = float(line[14:20])

  # Skip blank line
  fid.readline()

  # Get number density
  line = fid.readline()
  num_dens = float(line[22:34])

  # Get effective radius
  line = fid.readline()
  radius_eff = float(line[23:34])
  if l_humidity:
    line = fid.readline()
    radius_eff_dry = float(line[39:51])

  # Get volume fraction
  line = fid.readline()
  vol_frac = float(line[22:34])
  if l_humidity:
    line = fid.readline()
    vol_frac_dry = float(line[39:51])

  # Skip blank line and header
  for i in arange(3): fid.readline()

  # Read scattering properties
  band = zeros(0, dtype='i4')
  k_abs = zeros(0)
  k_scat = zeros(0)
  g_asym = zeros(0)
  n_newline = 0
  while True:
    line = fid.readline()

    # Check if end of data has been reached
    if n_newline > 1:
      # End of data
      break
    elif line == '\n':
      # Line is blank
      n_newline = n_newline + 1
    else:
      # Add data at current band
      band = append(band, int(line[:5]))
      k_abs = append(k_abs, float(line[9:25]))
      k_scat = append(k_scat, float(line[29:45]))
      g_asym = append(g_asym, float(line[49:65]))
      n_newline = 0

  # Convert absorption and scattering from m-1 to m2/kg
  if l_humidity:
    k_abs = k_abs/(vol_frac_dry*aer_density_component(aer_cmp))
    k_scat = k_scat/(vol_frac_dry*aer_density_component(aer_cmp))
  else:
    k_abs = k_abs/(vol_frac*aer_density_component(aer_cmp))
    k_scat = k_scat/(vol_frac*aer_density_component(aer_cmp))

  # All values of k_abs and k_scat should be >= 0, but check for negative
  # values just in case
  k_abs[k_abs < 0.0] = 0.0
  k_scat[k_scat < 0.0] = 0.0

  if l_humidity:
    return (aer_cmp, humidity, radius_eff, radius_eff_dry,
        band, k_abs, k_scat, g_asym)
  else:
    return aer_cmp, radius_eff, band, k_abs, k_scat, g_asym


def read_scatter_file(scatter_file, n_band, l_humidity, n_humidity=1):

  fid = open(scatter_file, 'r')

  if l_humidity:
    humidity = zeros(n_humidity)
    radius_eff = zeros(n_humidity)
    k_abs = zeros([n_humidity, n_band])
    k_scat = zeros([n_humidity, n_band])
    g_asym = zeros([n_humidity, n_band])
    for i_humid in arange(n_humidity):
      (aer_cmp, humidity[i_humid], radius_eff[i_humid], radius_eff_dry,
          band, k_abs[i_humid,:], k_scat[i_humid,:], g_asym[i_humid,:]) = (
          read_scatter_data(fid, l_humidity))

    fid.close()
    return (aer_cmp, humidity, radius_eff, radius_eff_dry,
        band, k_abs, k_scat, g_asym)

  else:
    aer_cmp, radius_eff_dry, band, k_abs, k_scat, g_asym = (
        read_scatter_data(fid, l_humidity))

    fid.close()
    return aer_cmp, radius_eff_dry, band, k_abs, k_scat, g_asym

################################################################################
# Begin main code
################################################################################

# Check if spectral file is for short-wave or long-wave
if "_lw_" in spectral_file_name:
  l_diag = True
else:
  l_diag = False

# Find number of radii
n_radius = zeros(len(aer_type), dtype='i4')
for i_aer in arange(len(aer_type)):
  radius = get_radius(scatter_data_dir, aer_type[i_aer])
  n_radius[i_aer] = len(radius)
n_radius_max = max(n_radius)

# Find number of humidities
l_humidity = zeros(len(aer_type), dtype='bool_')
n_humidity = zeros(len(aer_type), dtype='i4')
for i_aer in arange(len(aer_type)):
  radius = get_radius(scatter_data_dir, aer_type[i_aer])
  aer_file = os.path.join(scatter_data_dir,
      aer_type[i_aer] + '_' + radius[0] + '.avg')
  l_humidity[i_aer], n_humidity[i_aer] = humidity_dependence(aer_file)
n_humidity_max = max(n_humidity)

# Find number of bands
radius = get_radius(scatter_data_dir, aer_type[0])
aer_file = os.path.join(scatter_data_dir,
    aer_type[0] + '_' + radius[0] + '.avg')
n_band = get_n_band(aer_file, l_humidity[0])

# Read aerosol optical properties from .avg files
aer_cmp = zeros(len(aer_type), dtype='i4')
radius_eff_dry = zeros([len(aer_type), n_radius_max])
radius_eff = zeros([len(aer_type), n_radius_max, n_humidity_max])
band = zeros(n_band)
k_abs = zeros([len(aer_type), n_radius_max, n_humidity_max, n_band])
k_scat = zeros([len(aer_type), n_radius_max, n_humidity_max, n_band])
g_asym = zeros([len(aer_type), n_radius_max, n_humidity_max, n_band])
k_abs_550nm = zeros([len(aer_type), n_radius_max, n_humidity_max, 1])
k_scat_550nm = zeros([len(aer_type), n_radius_max, n_humidity_max, 1])
g_asym_550nm = zeros([len(aer_type), n_radius_max, n_humidity_max, 1])
for i_aer in arange(len(aer_type)):
  radius = get_radius(scatter_data_dir, aer_type[i_aer])

  for i_radius in arange(len(radius)):
    aer_file = os.path.join(scatter_data_dir,
        aer_type[i_aer] + '_' + radius[i_radius] + '.avg')
    if l_diag:
      aer_file_550nm = os.path.join(scatter_data_550nm_dir,
          aer_type[i_aer] + '_' + radius[i_radius] + '.avg')

    l_humidity[i_aer], n_humidity[i_aer] = humidity_dependence(aer_file)
    if l_humidity[i_aer]:
      (aer_cmp[i_aer], humidity, radius_eff[i_aer,i_radius],
          radius_eff_dry[i_aer,i_radius], band, k_abs[i_aer,i_radius,:,:],
          k_scat[i_aer,i_radius,:,:], g_asym[i_aer,i_radius,:,:]) = (
          read_scatter_file(aer_file, n_band,
          l_humidity[i_aer], n_humidity=n_humidity[i_aer]))
      if l_diag:
        (aer_cmp_tmp, humidity_tmp, radius_eff_tmp,
            radius_eff_dry_tmp, band_tmp, k_abs_550nm[i_aer,i_radius,:,:],
            k_scat_550nm[i_aer,i_radius,:,:],
            g_asym_550nm[i_aer,i_radius,:,:]) = (
            read_scatter_file(aer_file_550nm, 1, l_humidity[i_aer],
            n_humidity=n_humidity[i_aer]))
    else:
      (aer_cmp[i_aer], radius_eff_dry[i_aer,i_radius],
          band, k_abs[i_aer,i_radius,0,:], k_scat[i_aer,i_radius,0,:],
          g_asym[i_aer,i_radius,0,:]) = (
          read_scatter_file(aer_file, n_band, l_humidity[i_aer]))
      if l_diag:
        (aer_cmp_tmp, radius_eff_dry_tmp,
            band_tmp, k_abs_550nm[i_aer,i_radius,0,:],
            k_scat_550nm[i_aer,i_radius,0,:],
            g_asym_550nm[i_aer,i_radius,0,:]) = (
            read_scatter_file(aer_file_550nm, 1, l_humidity[i_aer]))

# Create nc file and add necessary dimensions
fout = Dataset(aer_nc, 'w', format='NETCDF3_CLASSIC')
fout.createDimension('component', len(unique(aer_cmp)))
fout.createDimension('radius', n_radius_max)
fout.createDimension('humidity', n_humidity_max)
fout.createDimension('band', n_band)

# Add dimension variables
component_nc = fout.createVariable('component', 'i4',
    dimensions=('component',))
radius_eff_nc = fout.createVariable('radius_eff', 'f8',
    dimensions=('component', 'radius', 'humidity'))
humidity_nc = fout.createVariable('humidity', 'f8',
    dimensions=('humidity',))
band_nc = fout.createVariable('band', 'i4',
    dimensions=('band',))
n_radius_nc = fout.createVariable('n_radius', 'i4',
    dimensions=('component',))
radius_eff_dry_nc = fout.createVariable('radius_eff_dry', 'f8',
    dimensions=('component', 'radius'))

# Add attributes to dimension variables
component_nc.title = 'aerosol component'
component_nc.long_name = 'aerosol component index as defined in rad_pcf'
radius_eff_dry_nc.title = 'effective dry radius'
radius_eff_dry_nc.long_name = 'effective radius of dry particle'
radius_eff_dry_nc.units = 'm'
radius_eff_nc.title = 'effective radius'
radius_eff_nc.long_name = 'effective radius of particle'
radius_eff_nc.units = 'm'
humidity_nc.title = 'humidities'
humidity_nc.long_name = 'humidities for which optical properties are tabulated'
band_nc.title = 'band'
band_nc.long_name = 'band'

# Add data variables
l_humidity_nc = fout.createVariable('l_humidity', 'i4',
    dimensions=('component',))
k_abs_nc = fout.createVariable('k_abs', 'f8',
    dimensions=('component', 'radius', 'humidity', 'band'))
k_scat_nc = fout.createVariable('k_scat', 'f8',
    dimensions=('component', 'radius', 'humidity', 'band'))
g_asym_nc = fout.createVariable('g_asym', 'f8',
    dimensions=('component', 'radius', 'humidity', 'band'))
density_nc = fout.createVariable('density', 'f8',
    dimensions=('component',))

# If short-wave file add variables for calculation of tau at 550 nm for
# diagnostics
if l_diag:
  k_abs_550nm_nc = fout.createVariable('k_abs_550nm', 'f8',
    dimensions=('component', 'radius', 'humidity'))
  k_scat_550nm_nc = fout.createVariable('k_scat_550nm', 'f8',
    dimensions=('component', 'radius', 'humidity'))
  g_asym_550nm_nc = fout.createVariable('g_asym_550nm', 'f8',
    dimensions=('component', 'radius', 'humidity'))

# Add attributes to data variables
l_humidity_nc.title = 'logical for humidity dependence'
l_humidity_nc.long_name = 'logical for inclusion of humidity dependnece'
k_abs_nc.title = 'absorption coefficient'
k_abs_nc.long_name = 'absorption coefficient'
k_abs_nc.units = 'm2/kg'
k_scat_nc.title = 'scattering coefficient'
k_scat_nc.long_name = 'scattering coefficient'
k_scat_nc.units = 'm2/kg'
g_asym_nc.title = 'asymmetry parameter'
g_asym_nc.long_name = 'asymmetry parameter'
density_nc.title = 'bulk density'
density_nc.long_name = 'bulk density of aerosol component'
if l_diag:
  k_abs_550nm_nc.title = 'absorption coefficient at 550 nm'
  k_abs_550nm_nc.long_name = 'absorption coefficient at 550 nm for diagnostics'
  k_abs_550nm_nc.units = 'm2/kg'
  k_scat_550nm_nc.title = 'scattering coefficient at 550 nm'
  k_scat_550nm_nc.long_name = 'scattering coefficient at 550 nm for diagnostics'
  k_scat_550nm_nc.units = 'm2/kg'
  g_asym_550nm_nc.title = 'asymmetry parameter at 550 nm'
  g_asym_550nm_nc.long_name = 'asymmetry parameter at 550 nm for diagnostics'

# Get indices of unique aerosol components
aer_cmp_unique, i_aer_cmp_unique = unique(aer_cmp, return_index=True)

# Set arrays common for all aerosol components
band_nc[:] = band
component_nc[:] = aer_cmp[i_aer_cmp_unique]
n_radius_nc[:] = n_radius[i_aer_cmp_unique]
l_humidity_nc[:] = l_humidity[i_aer_cmp_unique]
for i in arange(len(aer_cmp_unique)):
  density_nc[i] = aer_density_component(aer_cmp_unique[i])

# If humidity does not exist there are no aerosols with humidity dependence
if 'humidity' in locals():
  humidity_nc[:] = humidity # Assumes all aerosols have same humidity table
else:
  humidity_nc[:] = 0.0

# Add data for all aerosols
for i in arange(len(aer_cmp_unique)):
  i_aer = i_aer_cmp_unique[i]
  radius_eff_dry_nc[i,:] = (
      radius_eff_dry[i_aer,:])

  # Average aerosols with the same type index in aer_cmp
  i_aer_average = argwhere(aer_cmp_unique[i] == aer_cmp).flatten()

  # Check for humidity dependence as radius_eff, k_abs, k_scat and g_asym
  # depend on humidity
  if l_humidity[i_aer]:
    # Aerosol properties are humidity dependent
    radius_eff_nc[i,:,:] = radius_eff[i_aer_average,:,:]
    k_abs_nc[i,:,:,:] = mean(k_abs[i_aer_average,:,:,:], axis=0)
    k_scat_nc[i,:,:,:] = mean(k_scat[i_aer_average,:,:,:], axis=0)
    g_asym_nc[i,:,:,:] = mean(g_asym[i_aer_average,:,:,:], axis=0)
    if l_diag:
      k_abs_550nm_nc[i,:,:] = mean(k_abs_550nm[i_aer_average,:,:,0], axis=0)
      k_scat_550nm_nc[i,:,:] = mean(k_scat_550nm[i_aer_average,:,:,0], axis=0)
      g_asym_550nm_nc[i,:,:] = mean(g_asym_550nm[i_aer_average,:,:,0], axis=0)
  else:
    # No humidity dependence, use same values for all humidities
    for i_humid in arange(max(n_humidity)):
      radius_eff_nc[i,:,i_humid] = radius_eff_dry[i_aer]
      k_abs_nc[i,:,i_humid,:] = mean(k_abs[i_aer_average,:,0,:], axis=0)
      k_scat_nc[i,:,i_humid,:] = mean(k_scat[i_aer_average,:,0,:], axis=0)
      g_asym_nc[i,:,i_humid,:] = mean(g_asym[i_aer_average,:,0,:], axis=0)
      if l_diag:
        k_abs_550nm_nc[i,:,i_humid] = mean(k_abs_550nm[i_aer_average,:,0,0],
            axis=0)
        k_scat_550nm_nc[i,:,i_humid] = mean(k_scat_550nm[i_aer_average,:,0,0],
            axis=0)
        g_asym_550nm_nc[i,:,i_humid] = mean(g_asym_550nm[i_aer_average,:,0,0],
            axis=0)

# Add global attributes
fout.spectral_file = spectral_file_name
if l_diag:
  fout.spectral_region = 'lw'
else:
  fout.spectral_region = 'sw'

fout.close()

print('Successfully created ' + aer_nc)
print('All done')
