# Python script to generate Met Office P,T file for use when calculating
# absorption coefficients and k-coefficients
from pylab import *

# For hot Jupiters:
# T_min = 70.0 # Minimum temperature
# T_max = 3000.0 # Maximum temperature
# P_min = 1.0e-1 # Minimum pressure
# P_max = 1.0e+8 # Maximum pressure
# n_T = 20 # Number of temperature points
# n_P = 30 # Number of pressure points
# log_T = True
# file_header = ('File setting a standard set of {0} '.format(n_T*n_P) +
#     'combinations of pressure and temperature\n' +
#     'used to generate gaseous transmission data using correlated-k' +
#     'methods.'.format(n_T*n_P))

# For terrestrial planets
T_min = 100.0
T_max = 400.0
P_min = 1.0e+0
P_max = 1.0e+5
n_T = 13
n_P = int(round((log10(P_max) + log10(P_min))*10.0)) + 1
log_T = False
file_header = 'Write header here.'

# Generate P,T arrays
if log_T:
  T = logspace(log10(T_min), log10(T_max), n_T)
else:
  T = linspace(T_min, T_max, n_T)
P = logspace(log10(P_min), log10(P_max), n_P)

# Initialise writing of P,T file
f = open('pt{0}'.format(n_T*n_P), 'w')

# Write header
f.write(file_header + '\n')

# Write P,T table
f.write('*PTVAL')
for i in range(n_P):
  # Write pressure to file
  f.write('\n{0:13.7E}'.format(P[i]).replace('E', 'D'))
  for j in range(n_T):
    f.write(' {0:13.7E}'.format(T[j]).replace('E', 'D'))

f.write('\n*END')

# Close file
f.close()
