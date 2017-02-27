from pylab import *

def read_cld_data(spec_file, cld_type = 'drop'):

  if cld_type == 'drop':
    n_param = 16
    block_id = '*BLOCK: TYPE =   10'
  elif cld_type == 'ice':
    n_param = 14
    block_id = '*BLOCK: TYPE =   12'
  else:
    raise ValueError('Cloud type not recognised.')

  # Find drop/ice data
  fin = open(spec_file, 'r')
  while True:
    line = fin.readline()
    if line[:19] == block_id:
      block_version = int(line[47:51])
      break

  # Skip header
  for i in arange(4): fin.readline()

  # Get validity range of parametrization
  line = fin.readline()
  char_dim_min = float(line[40:52])
  char_dim_max = float(line[56:68])

  # Read data for each band
  param = zeros([n_band, n_param])
  fin.readline()
  for i in arange(n_band):
    line_long = ''
    while True:
      line = fin.readline()
      if line[:4] == 'Band' or line[:4] == '*END':
        break
      else:
        line_long = line_long + line
    param[i,:] = array(line_long.split())
  
  return param, char_dim_min, char_dim_max

def eval_drop_param(param, char_dim):
  k_ext = ((param[0] + char_dim*(param[1] + char_dim*param[2]))/
      (1.0 + char_dim*(param[3] + char_dim*(param[4] + char_dim*param[5]))))
  k_scat = (1.0 - (param[6] + char_dim*(param[7]+char_dim*param[8]))/
      (1.0 + char_dim*(param[9] + char_dim*param[10])))
  g_asym = ((param[11] + char_dim*(param[12] + char_dim*param[13]))/(1.0 +
      char_dim*(param[14] + char_dim*param[15])))

  return k_ext, k_scat, g_asym

def eval_ice_param(param, char_dim):
  x=param[3]/char_dim
  k_ext = (param[2]*x + param[1])*x + param[0]
  x = char_dim/param[8]
  k_scat = k_ext*(1.0 - (param[4] + x*(param[5] +
      x*(param[6] + x*param[7]))))
  x=char_dim/param[13]
  g_asym = param[9] + x*(param[10] + x*(param[11] + x*param[12]))

  return k_ext, k_scat, g_asym

def check_valid(char_dim, k_ext, k_scat, g_asym):
  param_ok = True
  if (any(k_ext < 0.0)):
    print('k_ext < 0.0')
    print('char_dim', char_dim)
    print('k_ext', k_ext)
    param_ok = False
  if (any(k_scat < 0.0)):
    print('k_scat < 0.0')
    print('char_dim', char_dim)
    print('k_scat', k_scat)
    param_ok = False
  if any(k_ext < k_scat):
    print('k_ext < k_scat')
    param_ok = False
  if any(g_asym > 1.0) or any(g_asym < 0.0):
    print('g_asym > 1.0 or g_asym < 0.0')
    print('char_dim', char_dim)
    print('g_asym', g_asym)
    param_ok = False

  return param_ok

# spec_file = '/home/damundse/Spectral_files/sp_lw_dsa_arcc/sp_lw_350_dsa_arcc'
spec_file = '/home/damundse/Spectral_files/sp_sw_dsa_ar/sp_sw_280_dsa_ar_trappist1'

n_point = 1000
plot_on = False

# Find number of spectral bands
fin = open(spec_file, 'r')
while True:
  line = fin.readline()
  if line[:19] == '*BLOCK: TYPE =    0':
    break
fin.readline()
line = fin.readline()
n_band = int(line[27:33])
fin.close()

drop_param, drop_char_dim_min, drop_char_dim_max = read_cld_data(spec_file,
    cld_type = 'drop')
ice_param, ice_char_dim_min, ice_char_dim_max = read_cld_data(spec_file,
    cld_type = 'ice')

drop_char_dim = logspace(log10(drop_char_dim_min),
    log10(drop_char_dim_max), n_point)
ice_char_dim = logspace(log10(ice_char_dim_min),
    log10(ice_char_dim_max), n_point)

all_ok = True
for i in arange(n_band):
  k_ext_drop, k_scat_drop, g_asym_drop = eval_drop_param(drop_param[i,:],
      drop_char_dim)
  k_ext_ice, k_scat_ice, g_asym_ice = eval_ice_param(ice_param[i,:],
      ice_char_dim)

  print('Band {:g}'.format(i+1))

  drop_ok = check_valid(drop_char_dim, k_ext_drop, k_scat_drop, g_asym_drop)
  ice_ok = check_valid(ice_char_dim, k_ext_ice, k_scat_ice, g_asym_ice)

  print('Drop: {}, Ice: {}'.format(drop_ok, ice_ok))

  if not drop_ok or not ice_ok:
    all_ok = False

  if plot_on:
    figure(1)
    loglog(drop_char_dim, k_ext_drop)
    figure(2)
    loglog(ice_char_dim, k_ext_ice)

if all_ok:
  print('All OK')
else:
  print('There are bad parameterisaitons')

if plot_on:
  figure(1)
  xlim([drop_char_dim_min, drop_char_dim_max])
  figure(2)
  xlim([ice_char_dim_min, ice_char_dim_max])
  show()
