import os

def read_refract_line(fin):

  line = fin.readline()

  # Check if end of data has been reached
  if line[3:8] == 'NSD=1':
    return -1, -1, -1
  
  # If end has not been reached proceed with reading
  else:
    wl = float(line[4:11])*1.0e-6
    re_refract = float(line[14:19])
    im_refract = float(line[20:29].replace('D', 'e'))
    
    if abs(wl - 9.9999e-4) < 1e-10:
      wl = 1.0e-3
    
    if im_refract < 1.0e-8:
      im_refract = 1.0e-8
    
    return wl, re_refract, im_refract

def write_refract_line(fout, wl, re_refract, im_refract):
  fout.write('{wl:16.5e}{re:16.5e}{im:16.5e}\n'.format(
      wl=wl, re=re_refract, im=im_refract))

def convert_mieqpc_socrates(file_mieqpc, file_socrates):

  fin = open(file_mieqpc, 'r')
  fout = open(file_socrates, 'w')

  # Write header
  fout.write('\n')
  fout.write('     Wavelength (m)  Real Part       Imaginary Part\n')
  fout.write('*BEGIN_DATA\n')

  # Read first line of data
  fin.readline()
  wl, re_refract, im_refract = read_refract_line(fin)

  # If wavelength is > 0.2 micron, extrapolate to 0.2 micron
  if wl > 2.0e-7:
    write_refract_line(fout, 2.0e-7, re_refract, im_refract)
  write_refract_line(fout, wl, re_refract, im_refract)

  # Loop to read over all wavelengths
  while True:
    wl, re_refract, im_refract = read_refract_line(fin)
    
    # Check if end of file has been reached
    if wl < 0.0:
      break

    # If not end of file, write data
    else:
      write_refract_line(fout, wl, re_refract, im_refract)

  # Set refractive index at 0.1 m to (1,0)
  write_refract_line(fout, 0.1, 1.0, 0.0)

  fout.write('*END')
  fin.close()
  fout.close()
  
path_mieqpc = '/Users/damundse/Software/GISS_Rad/gcm_table/datamie'
files_mieqpc = [
    'mieqpc.22antwater',
    'mieqpc.22organic',
    'mieqpc.22seasalt',
    'mieqpc.22sulfate',
    'mieqpc.h2so4_v1',
    'mieqpc.soot',
    'qmie25.sdust',
    'mieqpc.sinyukdust',
    ]

path_socrates = '/Users/damundse/Scatter_data'
files_socrates = [
    'refract_antwater',
    'refract_organic',
    'refract_seasalt',
    'refract_sulfate',
    'refract_h2so4',
    'refract_soot',
    'refract_sdust',
    'refract_sinyukdust',
    ]

for i in range(len(files_mieqpc)):
  convert_mieqpc_socrates(os.path.join(path_mieqpc, files_mieqpc[i]),
      os.path.join(path_socrates, files_socrates[i]))
