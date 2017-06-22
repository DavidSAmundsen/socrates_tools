from pylab import *
from scipy.interpolate import interp2d

file_in = 'CO2-H2_200_250_300_350.cia'
file_out = 'h2-co2.cia'
chem_sym = 'H2-CO2'

# file_in = 'CO2-CH4_200_250_300_350.cia'
# file_out = 'co2-ch4.cia'
# chem_sym = 'CO2-CH4'

# Define reference to data source
reference = r'Wordsworth+ GRL 2017'

# Load data
data = loadtxt(file_in)

# Convert CIA from [cm^-1/amagat^2] to [cm^5/molecule^2]
loschmidt_cnst = 2.6867774e+19
data[:,1:] = data[:,1:]/(loschmidt_cnst**2)

# Set pressure and wavenumber grids
t_file = array([200.0, 250.0, 300.0, 350.0])
t_cia = arange(min(t_file), max(t_file) + 10.0, 10.0)
nu_step = 1.0
nu_cia = arange(nu_step, max(data[:,0]) + nu_step, nu_step)

# Check that number of temperatures agrees with number of data columns
n_t_file = size(data, 1) - 1
if n_t_file != len(t_file):
  raise NameError('Temperatures do not agree.')

# If first wavenumber is 0.0, replace with small value
if all(data[0,:] == 0.0):
  data[0,:] = data[1,:]*1e-10

# Perform interpolation in log T, log CIA, and linear in wavenumber
fip = interp2d(log10(t_file), data[:,0], log10(data[:,1:]))
cia = 10.0**fip(log10(t_cia), nu_cia)

# Loop through temperatures and write CIA on HITRAN format
fout = open(file_out, 'w')
for i in arange(len(t_cia)):
  header = ('{chem_sym:20s}'.format(chem_sym=chem_sym) +
            '{nu_min:10.3f}'.format(nu_min=nu_cia[0]) +
            '{nu_max:10.3f}'.format(nu_max=nu_cia[-1]) +
            '{n_nu:7g}'.format(n_nu=len(nu_cia)) +
            '{tmp:7.1f}'.format(tmp=t_cia[i]) +
            ' '*22 +
            '{ref:21s}'.format(ref=reference) +
            '\n')
  fout.write(header)
  for j in arange(len(nu_cia)):
    line = '{nu:10.3f} {cia:10.3E}\n'.format(nu=nu_cia[j], cia=cia[j,i])
    fout.write(line)
fout.close()
