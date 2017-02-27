from pylab import *

fit_file = '/home/damundse/k-coeff/sw_dsa_ar/cld/sun/fit_sw_drop5_29'
mon_file = '/home/damundse/k-coeff/sw_dsa_ar/cld/sun/mon_sw_drop5_29'
n_param = 16
n_band = 29
n_point = 52
n_fit = 3
r_eff_min = 1.7e-06
r_eff_max = 5.0e-05

colors = ['b', 'g', 'r', 'c', 'm', 'y',
          'b', 'g', 'r', 'c', 'm', 'y',
          'b', 'g', 'r', 'c', 'm', 'y',
          'b', 'g', 'r', 'c', 'm', 'y',
          'b', 'g', 'r', 'c', 'm', 'y',
          'b', 'g', 'r', 'c', 'm', 'y',
          'b', 'g', 'r', 'c', 'm', 'y']

# Read final fits from fit file
i_band = 0
i_param = 0
read_fit = False
cld_param = zeros([n_band, n_param])
for line in open(fit_file, 'r'):
  if not read_fit and line[0:4] == 'Band':
    read_fit = True
    continue
  if read_fit:
    cld_param_str = line[4:].replace(' '*4, ',').replace('\n', '').split(',')
    cld_param_lst = [float(str) for str in cld_param_str]
    cld_param[i_band, i_param:i_param+4] = array(cld_param_lst)
    i_param = i_param + 4
    if i_param == n_param:
      read_fit = False
      i_band = i_band + 1
      i_param = 0

# Read goodness of fit from the monitoring file
i_band = 0
i_point = 0
i_fit = 0
read_fit = False
skip_line = False
cld_fit = zeros([n_band, n_fit, n_point, 3])
for line in open(mon_file, 'r'):
  if skip_line:
    skip_line = False
    continue
  if not read_fit and line[4:15] == 'Scaled Size':
    read_fit = True
    skip_line = True
    continue
  if read_fit:
    cld_fit_str = line[3:].replace(' '*3, ',').replace('\n', '').split(',')
    cld_fit_lst = [float(str) for str in cld_fit_str]
    cld_fit[i_band, i_fit, i_point, :] = array(cld_fit_lst)
    i_point = i_point + 1
    if i_point == n_point:
      read_fit = False
      i_point = 0
      i_fit = i_fit + 1
      if i_fit == n_fit:
        i_fit = 0
        i_band = i_band + 1

# Plot actual data and fits for each band
for i_fit in arange(n_fit):
  figure(i_fit)
  for i_band in arange(n_band):
    semilogx(cld_fit[i_band, i_fit, :, 0], cld_fit[i_band, i_fit, :, 1], 
        '.' + colors[i_band])
    semilogx(cld_fit[i_band, i_fit, :, 0], cld_fit[i_band, i_fit, :, 2], 
        '-' + colors[i_band])

# Plot k_ext, k_scat and g at high resolution in r_eff
r_eff = logspace(log10(r_eff_min), log10(r_eff_max), 1000)
for i_band in arange(n_band):
  k_ext = (cld_param[i_band,0]+r_eff*(cld_param[i_band,1]+r_eff
      *cld_param[i_band,2]))/(1.0e+00+r_eff*(cld_param[i_band,3]+r_eff
      *(cld_param[i_band,4]+r_eff*cld_param[i_band,5])))
  k_scat = k_ext*(1.0e+00-(cld_param[i_band,6]+r_eff*(cld_param[i_band,7]
      +r_eff*cld_param[i_band,8]))
      /(1.0e+00+r_eff*(cld_param[i_band,9]+r_eff*cld_param[i_band,10])))
  asym = (cld_param[i_band,11]+r_eff*(cld_param[i_band,12]+r_eff
      *cld_param[i_band,13]))/(1.0e+00+r_eff
      *(cld_param[i_band,14]+r_eff*cld_param[i_band,15]))
  figure(n_fit)
  semilogx(r_eff, k_ext)
  figure(n_fit+1)
  semilogx(r_eff, k_scat)
  figure(n_fit+2)
  semilogx(r_eff, asym)

show()
