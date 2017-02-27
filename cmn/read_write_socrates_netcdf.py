# Functions to create and write SOCRATES NetCDF files for use with Cl_run_cdf.

from netCDF4 import *
import os

# Writes some quantity x as a function of pressure p_lev
def write_latlon(path, x, prefix, suffix, name, unit, title,
    p_lev=None, basis=None):

  if p_lev is not None and basis is not None:
    print('Error: Both p_lev and basis brovided.')
    return
  
  path_file = os.path.join(path, prefix + '.' + suffix)
  
# Create the output netCDF file and get output file id
  fid = Dataset(path_file, mode='w', format='NETCDF3_CLASSIC')
  
# Add dimensions to the file
  fid.createDimension('lon', 1)
  fid.createDimension('lat', 1)
  if p_lev is not None:
    fid.createDimension('plev', len(p_lev))
  elif basis is not None:
    fid.createDimension('basis', len(basis))

# Create variables and attributes
  lon_nc = fid.createVariable('lon', 'f8', dimensions=('lon'))
  lon_nc.units = 'degree'
  lon_nc.title = 'Latitude'
  lon_nc[:] = 0.0
  
  lat_nc = fid.createVariable('lat', 'f8', dimensions=('lat'))
  lat_nc.units = 'degree'
  lat_nc.title = 'Latitude'
  lat_nc[:] = 0.0
  
  if p_lev is None and basis is None:
    x_nc = fid.createVariable(suffix, 'f8', dimensions=('lat', 'lon'))
  
  elif p_lev is not None:
    plev_nc = fid.createVariable('plev', 'f8', dimensions=('plev'))
    plev_nc.units = 'Pa'
    plev_nc.title = 'Pressure'
    
    x_nc = fid.createVariable(suffix, 'f8', dimensions=('plev', 'lat', 'lon'))
  elif basis is not None:
    basis_nc = fid.createVariable('basis', 'i4', dimensions=('basis'))
    basis_nc.units = 'None'
    basis_nc.title = "Basis function"
  
    x_nc = fid.createVariable(suffix, 'f8', dimensions=('basis', 'lat', 'lon'))
    
  x_nc.units = unit
  x_nc.title = title
  
# Write arrays to file
  if p_lev is not None:
    plev_nc[:] = p_lev
  elif basis is not None:
    basis_nc[:] = basis
  x_nc[:] = x
  
  fid.close()
  

# Read some quantity from a SOCRATES netCDF file
def read_latlon(path, prefix, suffix, name):

  path_file = os.path.join(path, prefix + '.' + suffix)
  
# Open netCDF file
  fid = Dataset(path_file, mode='r')
  
# Read pressure and requested quantity
  plev = fid.variables['plev'][:]
  x = fid.variables[name][:,0,0]
  
  return plev, x