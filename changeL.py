import os
import sys

LMIN = 6.75
LMAX = 8.0

# -----------------
# Config file:
config_file = str(sys.argv[1])
counter = float(str(sys.argv[2]))
loops = float(str(sys.argv[3]))

EL = (LMAX-LMIN)*counter/loops + LMIN
text = ''

with open(config_file, 'r') as f:
	for line in f:
		if line[0] == 'L':
			text += 'L ' + str(EL) +'\n'
			print 'L = ' + str(EL) + ' fm'
		else:
			text += line

with open(config_file+'.tmp', 'w') as f:
	f.write(text)

os.rename(config_file+'.tmp',config_file)