# Removes negative entries from HITRAN CIA data.

file_in = r'/scr2/socrates/continua/n2-n2_2011.cia'
file_out = r'/scr2/socrates/continua/n2-n2_2011_clean.cia'

# Count negative entries
fid_in = open(file_in, 'r')
line_in = fid_in.readline()
n_entry = 1
n_neg = [0]
while True:
  n_nu = int(line_in[40:47])
  for i in range(n_nu):
    line_in = fid_in.readline()
    nu = float(line_in[:10])
    kabs = float(line_in[11:21])
    if kabs < 0.0:
      n_neg[n_entry-1] = n_neg[n_entry-1] + 1
  line_in = fid_in.readline()
  if line_in == '':
    break
  else:
    n_entry = n_entry + 1
    n_neg.append(0)
fid_in.close()

print(n_neg)

# Remove negative entries
fid_in = open(file_in, 'r')
fid_out = open(file_out, 'w')
for i_entry in range(n_entry):
  line_in = fid_in.readline()
  n_nu = int(line_in[40:47])
  line_out = (line_in[:40] +
      '{n_nu:7g}'.format(n_nu = n_nu - n_neg[i_entry]) +
      line_in[47:])
  fid_out.write(line_out)
  for i in range(n_nu):
    line_in = fid_in.readline()
    nu = float(line_in[:10])
    kabs = float(line_in[11:21])
    if kabs >= 0.0:
      line_out = '{nu:10.4f} {kabs:10.3E}\n'.format(nu = nu, kabs = kabs)
      fid_out.write(line_out)

fid_in.close()
fid_out.close()
