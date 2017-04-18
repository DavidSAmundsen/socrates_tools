from pylab import *

spec_file = sys.argv[1]

# Get number of spectral bands
fin = open(spec_file, 'r')
while True:
  line = fin.readline()
  if line[:27] == 'Number of spectral bands = ':
    n_band = int(line[27:])
    break
fin.close()

# Get number of absorbing gases
fin = open(spec_file, 'r')
while True:
  line = fin.readline()
  if line[:36] == 'Total number of gaseous absorbers = ':
    n_gas = int(line[36:])
    break
fin.close()

# Create array for number of k-terms
n_k_terms = zeros([n_band, n_gas], dtype='int')

# Find Block 5
fin = open(spec_file, 'r')
while True:
  line = fin.readline()
  if line[:19] == '*BLOCK: TYPE =    5':
    break

# Skip Block 5 header
for i in arange(4):
  fin.readline()

# Read data for all bands and all gases
while True:
  line = fin.readline()
  
  # Check if end of Block 5 has been reached
  if line[:4] == '*END':
    break
  
  try:
    i_band = int(line[:5])
  except:
    continue
  i_gas = int(line[6:17])
  n_k_terms[i_band-1, i_gas-1] = int(line[18:29])
  i_scaling = int(line[42:53])

  if i_scaling == 2:
    lines_to_skip = 2*n_k_terms[i_band-1, i_gas-1] + 1
  else:
    lines_to_skip = n_k_terms[i_band-1, i_gas-1] + 1

  # Skip k-terms
  for i in arange(lines_to_skip):
    fin.readline()

# Print results
for i_gas in arange(n_gas):
  print('Gas ', i_gas + 1, ' k-terms:', int(sum(n_k_terms[:, i_gas])))
print('Total number of k-terms: ', int(sum(sum(n_k_terms))))

print(n_k_terms)

# Calculate number of monochromatic calculations required if using random
# overlap without resorting and rebinning
n_mono_ro = 0
for i_band in arange(n_band):
  n_mono_ro = n_mono_ro + prod(n_k_terms[i_band, n_k_terms[i_band,:] > 0])

print('Number of monchromatic calculations with RO', n_mono_ro)

fin.close()
