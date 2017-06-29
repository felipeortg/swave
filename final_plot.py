import os
import sys
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import scipy.optimize as opt
import datetime as dt
import numpy as np
import re

# -----------------
# Config file
config_file = str(sys.argv[1])
sd = False

with open(config_file, 'r') as f:
	for line in f:
		if line[0:6] == 'folder':
			folder = line[7:-1]
		if line[0] == 'L':
			L_lattice_size = float(line[2:]) # fermis
		if line[0] == 'l':
			l = int(line[2:])
		if line[0] == 'm':
			m = int(line[2:])
		if line[0:2] == 'dx':
			dx = int(line[2:])
		if line[0:2] == 'dy':
			dy = int(line[2:])
		if line[0:2] == 'dz':
			dz = int(line[2:])
		if line[0:5] == 'x2max':
			x2max = float(line[6:])
		if line[0:5] == 'x2min':
			x2min = float(line[6:])
		if line[0:6] == 'points':
			points = int(line[7:])
		if line[0:2] == 'g0':
			g0 = float(line[2:])
		if line[0:2] == 'M0':
			m0 = float(line[3:])

d = [dx,dy,dz]

# Define constants
I = complex(0.0,1.0)
hbar = 197.3269718 # MeV fm when speed of ligth is 1

# Define variables of the C function
# S-wave
# d= [0.0,0.0,0]
# L_lattice_size = 10.0 # fermis
# L = L_lattice_size/hbar # MeV^-1
# lab_moment = (2*np.pi/L)*np.array(d)
# lab_moment2 = np.dot(lab_moment,lab_moment)
m1 = 146.0 # MeV
m2 = 146.0 # MeV
# m0 = 770.0 # MeV, Resonance mass
# g0 = 3.0 # coupling

eles = [6.0,6.5,7.0,7.5,8.0,10.0]

# Matplotlib param
plt.rc('text', usetex=True)
plt.rc('font', size=12)
plt.rc('font', family='serif')
plt.rc('axes.formatter', useoffset = False)

specfolder = '../'+folder+'Spectrum_' + str(dt.datetime.now())[:-16] +'/'

fig = plt.figure()

for spectrum in os.listdir(specfolder):
	if spectrum[0] == '.' or spectrum[-1] == 'f':
		continue

	ele = float(re.findall('\d+.\d+',spectrum)[0])

	cm_energy_spec = np.load(specfolder+spectrum)

	plt.plot(ele*np.ones(len(cm_energy_spec)), cm_energy_spec,'o',color='r')

for nn in range(6):
	free_part_en = lambda Llat: 2*np.sqrt(nn*((2*np.pi*hbar)/Llat)**2 + m1**2)
	Lspace = np.linspace(5.5, 8, num=100)
	plt.plot(Lspace,free_part_en(Lspace),'--',color='k')

plt.plot([5.5,8],[m0,m0],'--',color='0.5')

plt.title(r'Spectrum with BW ($m_0$ = '+str(m0)+' MeV, $g_0$ = '+ str(g0) +')')

plt.xlabel(r'$L/$fm')

plt.ylabel(r'$E/$MeV')

figname = specfolder + 'Spectrum BW (m_0 = '+str(m0)+' MeV, g_0 = '+ str(g0) +')--' + str(dt.datetime.now())[:-16]+'.pdf'

fig.savefig(figname)

plt.show()
