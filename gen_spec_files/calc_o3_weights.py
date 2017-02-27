from pylab import *

def get_flux(spec_file, i_band):
  fin = open(spec_file, 'r')

  # Locate block 2 with solar spectrum
  while True:
    line = fin.readline()
    if line[:19] == '*BLOCK: TYPE =    2':
      break

  # Skip header
  for i in arange(2): fin.readline()

  # Go to band
  for i in arange(i_band): line = fin.readline()

  # Get flux in band
  flux = float(line[13:28])
  return flux

# 6 and 12 band spectral files are input arguments
if (len(sys.argv) <= 2):
  raise ValueError('Both 6 and 12 band spectral files must be provided')
spec_file_6 = sys.argv[1]
spec_file_12 = sys.argv[2]

# Array to keep weights
weight = zeros(8)

# Get weights for bands 1 to 6 in 12 band file, corresponding to band
# 1 in 6 band file
for i_band in arange(1, 7):
  weight[i_band-1] = get_flux(spec_file_12, i_band)/get_flux(spec_file_6, 1)

# Get weights for bands 7 and 8 in 12 band file, corresponding to band
# 2 in 6 band file
for i_band in arange(7, 9):
  weight[i_band-1] = get_flux(spec_file_12, i_band)/get_flux(spec_file_6, 2)

# Construct weights on required format
weight_frmt = zeros(8)
for i_band in arange(8):
  print(i_band+1, '{:.9E}'.format(weight[i_band]))