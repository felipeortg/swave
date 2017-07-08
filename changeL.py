import os
import sys

# -----------------
# Config file:
config_file = str(sys.argv[1])
counter = float(str(sys.argv[2]))
loops = float(str(sys.argv[3]))

with open(config_file, 'r') as f:
	for line in f:
		if line[0:5] == 'LMIN ':
			LMIN = float(line[5:])
		if line[0:5] == 'LMAX ':
			LMAX = float(line[5:])

EL = (LMAX-LMIN)*counter/(loops-1) + LMIN
text = ''

with open(config_file, 'r') as f:
	for line in f:
		if line[0:2] == 'L ':
			text += 'L ' + str(EL) +'\n'
			print 'L = ' + str(EL) + ' fm'
		else:
			text += line

with open(config_file+'.tmp', 'w') as f:
	f.write(text)

os.rename(config_file+'.tmp',config_file)
