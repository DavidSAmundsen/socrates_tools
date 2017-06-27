# Python script to construct NetCDF input files for use with Cl_run_cdf.

from pylab import *
import os

from read_write_socrates_netcdf import write_latlon
from gas_list import gas_name

################################################################################
# User input
################################################################################

# Charnay et al. Archean Case A/B/C
# run_path = '/home/damundse/socrates_runs/terrestrial/'
# from modele_writer_pt_profiles import pt_profile, defines_pressure
# from mass_mixing_ratio_archean import mass_mixing_ratio
# flag_mmr = 'Charnay_et_al_Case_C'
# run_name = 'ar_case_c'
# profile_name = ('/home/damundse/ModelE_Support/rad_data/' +
#     'Archean_midlat_summer/charnay_case_b_earth_tbl.dat')
# star_irr_toa = 1381.527
# star_zenith = 0.0
# albedo_lw = 0.0
# albedo_sw = 0.1
# basis = array([0])
# surf_temp = None
# absorbers = [
#     'H2O',
#     'CO2',
#     'CH4',
#     ]
# # mix_ratio_override = {
# #     'H2O': 4e-2,
# #     'CO2': 1e-2,
# #     'CH4': 0.0057
#     }
# p_surf = 1.0e+4

# Mars
# run_path = '/home/damundse/socrates_runs/terrestrial/'
# from modele_writer_pt_profiles import pt_profile, defines_pressure
# from mass_mixing_ratio_mars import mass_mixing_ratio
# flag_mmr = 'mars'
# run_name = 'mars'
# profile_name = ('/home/damundse/ModelE_Support/rad_data/' +
#     'Archean_midlat_summer/charnay_case_b_earth_tbl.dat')
# star_irr_toa = 1381.527
# star_zenith = 0.0
# albedo_lw = 0.0
# albedo_sw = 0.5
# basis = array([0])
# surf_temp = None
# absorbers = [
#     'H2O',
#     'CO2',
#     ]
# mix_ratio_override = {
# #     'H2O': 4e-2,
#     }
# p_surf = 1.0e+6

# Vlad Case A/B/C
run_path = '/home/damundse/socrates_runs/terrestrial/'
from modele_writer_pt_profiles import pt_profile, defines_pressure
from mass_mixing_ratio_archean import mass_mixing_ratio
flag_mmr = 'Vlad_Case_C'
run_name = 'vlad_case_c'
profile_name = ('/home/damundse/ModelE_Support/rad_data/' +
    'Archean_midlat_summer/charnay_case_b_earth_tbl.dat')
star_irr_toa = 1381.527
star_zenith = 0.0
albedo_lw = 0.0
albedo_sw = 0.0
basis = array([0])
surf_temp = None
absorbers = [
    'H2O',
    'CO2',
    'CH4',
    'N2O',
    ]

# Paleo Mars
# run_path = '/home/damundse/socrates_runs/terrestrial/'
# from modele_writer_pt_profiles import pt_profile, defines_pressure
# from mass_mixing_ratio_paleo_mars import mass_mixing_ratio
# flag_mmr = 'setup_4'
# run_name = 'pmars_setup4'
# profile_name = ('/home/damundse/ModelE_Support/rad_data/' +
#     'Archean_midlat_summer/charnay_case_b_earth_tbl.dat')
# star_irr_toa = 1381.527
# star_zenith = 0.0
# albedo_lw = 0.0
# albedo_sw = 0.5
# basis = array([0])
# surf_temp = None
# absorbers = [
#     'H2O',
#     'CO2',
#     'CH4',
#     'N2',
#     'H2',
#     ]
# mix_ratio_override = {
# #    'H2O': 0.0,
# #    'CO2': 0.0, # 1e-6*44./28.,
# #    'CH4': 1.4e-2*16./28.,
#     }
# p_surf = 1.0e+6

# Humid atmosphere
# run_path = '/home/damundse/socrates_runs/terrestrial/'
# from modele_writer_pt_profiles import pt_profile, defines_pressure
# from mass_mixing_ratio_archean import mass_mixing_ratio
# flag_mmr = 'Vlad_Case_B'
# run_name = 'humid'
# profile_name = ('/home/damundse/ModelE_Support/rad_data/' +
#     'Archean_midlat_summer/charnay_case_b_earth_tbl.dat')
# star_irr_toa = 1.2*1381.527
# star_zenith = 0.0
# albedo_lw = 0.0
# albedo_sw = 0.0
# basis = array([0])
# surf_temp = None
# absorbers = [
#     'H2O',
#     'CO2',
# #     'AIR',
#     ]
# mix_ratio_override = {
# #     'H2O': 0.1,
#     'CO2': 1e-6,
# #     'AIR': 0.96,
#     }
    
# US Standard Atmosphere 1976
# from modele_writer_pt_profiles import pt_profile, defines_pressure
# from mass_mixing_ratio_standard_earth import mass_mixing_ratio
# flag_mmr = 'tropical'
# run_name = 'tropical'
# profile_name = ('/Users/damundse/ModelE_Support/rad_data/' +
#     'Present_day_Earth/tropical.dat')
# star_irr_toa = 1381.527
# star_zenith = 30.0
# albedo_lw = 0.0
# albedo_sw = 0.0
# basis = array([0])
# surf_temp = None
# absorbers = [
#     'H2O',
#     'O3',
#     'CO2',
#     'N2O',
#     'CH4',
#     ]
    
# Midlatitude summer atmosphere
# from modele_writer_pt_profiles import pt_profile, defines_pressure
# from mass_mixing_ratio_standard_earth import mass_mixing_ratio
# run_path = '/home/damundse/socrates_runs/terrestrial/'
# flag_mmr = 'midlat_summer'
# run_name = 'midlat_summer'
# profile_name = ('/home/damundse/ModelE_Support/rad_data/' +
#     'Present_day_Earth/midlat_summer.dat')
# star_irr_toa = 1381.527
# star_zenith = 0.0
# albedo_lw = 0.0
# albedo_sw = 0.5
# basis = array([0])
# surf_temp = None
# absorbers = [
#     'H2O',
#     'O3',
#     'CO2',
#     'N2O',
#     'CH4',
#     'O2',
#     ]
# mix_ratio_override = {
#     'H2O': 0.01,
#     }
    
# Set PT profile manually
# p_lev_min = 1e+0
# p_lev_max = 1.013e5
# n_lay = 50
# include_ghost_lay = False
# defines_pressure = False
# def pt_profile(P, profile_name):
#   return 285.0
t_enhance = 1.06

################################################################################
# Get P-T profile
################################################################################

# Define flag_mmr for mass mixing routine if not defined
if 'flag_mmr' not in vars():
  flag_mmr = None

# Define mix_ratio_override if not defined
if 'mix_ratio_override' not in vars():
  mix_ratio_override = {
      }

# Print flag_mmr if provided
if flag_mmr:
  print('Using mass mixing ratio flag ' + str(flag_mmr))

# Crate directory if it does not exist
run_path = os.path.join(run_path, run_name)
if not os.path.isdir(run_path):
  os.mkdir(run_path)

# If PT profile function does not define the pressure
if (not defines_pressure):

# Construct logarithmic pressure arrays
  p_lev = logspace(log10(p_lev_min), log10(p_lev_max), n_lay+1)
  p_lay = (p_lev[1:] + p_lev[:-1])*0.5

# Add ghost layer if requested
  if include_ghost_lay:
    p_lev = insert(p_lev, 0, 0.0)
    p_lay = insert(p_lay, 0, p_lev[1]*0.5)
    i_p_start = 1
  else:
    i_p_start = 0

# Get temperatures from selected P,T profile
  t_lev = zeros(size(p_lev))
  t_lay = zeros(size(p_lay))
  for i in arange(i_p_start, len(p_lay)):
    t_lev[i] = pt_profile(p_lev[i], profile_name)
    t_lay[i] = pt_profile(p_lay[i], profile_name)
  t_lev[-1] = pt_profile(p_lev[-1], profile_name)

# Get temperature of ghost layer if requested
  if include_ghost_lay:
    t_lev[0] = 0.0
    t_lay[0] = pt_profile(p_lay[0], profile_name)

# PT profile function defines pressure
else:
  p_lev, t_lev, p_lay, t_lay = pt_profile(profile_name)

# Set surface temperature to temperature of lowest level if requested
if surf_temp is None:
  surf_temp = t_lev[-1]

# Enhance temperatures if requested
if 't_enhance' in locals():
  t_lev = t_lev*t_enhance
  t_lay = t_lay*t_enhance

################################################################################
# Get mass mixing ratios
################################################################################

mmr = dict()
for absorber in absorbers:
  if absorber in mix_ratio_override.keys():
    mmr[absorber] = mix_ratio_override[absorber]*ones(len(p_lay))
  else:
    mmr[absorber] = mass_mixing_ratio(absorber, p_lay, t_lay, flag = flag_mmr,
        return_dry = True)

################################################################################
# Extend to higher pressures if required
################################################################################

if 'p_surf' in locals():
  if p_surf > p_lev[-1]:
    p_enhance = p_surf/p_lev[-1]
    p_lev_new = p_lev*p_enhance
    i_add = argwhere(p_lev > p_lev_new[1])[0,0]-1
    p_lev = append(p_lev[:i_add], p_lev_new[1:])
    p_lay = (p_lev[1:] + p_lev[:-1])*0.5
    t_lev = append(t_lev[0]*ones(i_add), t_lev[1:])
    t_lay = append(t_lev[0]*ones(i_add), t_lay[1:])
    # Restrict temperature at model top
    if any(t_lev[:i_add] > t_lev[i_add]):
      t_lay[:i_add] = t_lev[i_add]*ones(i_add)
      t_lev[:i_add] = t_lev[i_add]*ones(i_add)
    for absorber in absorbers:
      mmr[absorber] = append(mmr[absorber][0]*ones(i_add), mmr[absorber][1:])
  elif p_surf < p_lev[-1]:
    p_lay = p_lay*p_surf/p_lev[-1]
    p_lev = p_lev*p_surf/p_lev[-1]

################################################################################
# Create NetCDF files
################################################################################

# Temperatures
write_latlon(run_path, t_lev, run_name, 'tl', 'tl', 'K',
    'Temperature on levels', p_lev=p_lev)
write_latlon(run_path, t_lay, run_name, 't', 't', 'K',
    'Temperature', p_lev=p_lay)

# Mass mixing ratios
for absorber in absorbers:
  if (absorber == 'H2O'):
    write_latlon(run_path, mmr[absorber], run_name, 'q',
        'q', 'None', gas_name[absorber] + ' MMR', p_lev=p_lay)
  else:
    write_latlon(run_path, mmr[absorber], run_name, absorber.lower(),
        absorber.lower(), 'None', gas_name[absorber] + ' MMR', p_lev=p_lay)

# Write albedos
write_latlon(run_path, albedo_lw, run_name, 'surf_lw', 'alb', 'None',
    'Albedo weights', basis=basis)
write_latlon(run_path, albedo_sw, run_name, 'surf_sw', 'alb', 'None',
    'Albedo weights', basis=basis)

# Irradiation at top of atmosphere
write_latlon(run_path, star_irr_toa, run_name, 'stoa', 'stoa', 'W m-2',
    'Star irradiance')

# Star zenith angle
write_latlon(run_path, star_zenith, run_name, 'szen', 'szen', 'degrees',
    'Star zenith angle', basis=basis)

# Surface temperature
write_latlon(run_path, surf_temp, run_name, 'tstar', 'tstar', 'K',
    'Surface temperature', p_lev=array([p_lev[-1]]))
