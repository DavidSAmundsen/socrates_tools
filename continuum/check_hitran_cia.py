# Checks if HITRAN CIA data contains negative entries.
import sys

file = sys.argv[1]

fid = open(file, 'r')

found_neg = False
line = fid.readline()
n_nu = int(line[40:47])
while True:
  n_nu = int(line[40:47])
  for i in range(n_nu):
    line = fid.readline()
    nu = float(line[:10])
    kabs = float(line[11:21])
    if kabs < 0.0:
      found_neg = True
      break
  line = fid.readline()
  if found_neg:
    break
  elif line == '':
    break

if found_neg:
  print('Negative values found.')
else:
  print('No negative values found. All OK.')