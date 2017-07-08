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
		if line[0:7] == 'folder ':
			folder = line[7:-1]
		if line[0:3] == 'dx ':
			dx = int(line[3:])
		if line[0:3] == 'dy ':
			dy = int(line[3:])
		if line[0:3] == 'dz ':
			dz = int(line[3:])
		if line[0:3] == 'g0 ':
			g0 = float(line[3:])
		if line[0:3] == 'M0 ':
			m0 = float(line[3:])
		if line[0:3] == 'sd ':
			sd = True
			speci_dat = str(line[3:-1])
		if line[0:5] == 'LMIN ':
			elemin = float(line[5:])
		if line[0:5] == 'LMAX ':
			elemax = float(line[5:])

d = [dx,dy,dz]

# Define constants
I = complex(0.0,1.0)
hbar = 197.3269718 # MeV fm when speed of ligth is 1

# Define variables of the C function

m1 = 146.0 # MeV
m2 = 146.0 # MeV


# Matplotlib param
plt.rc('text', usetex=True)
plt.rc('font', size=12)
plt.rc('font', family='serif')
plt.rc('axes.formatter', useoffset = False)

if sd == False:
	specfolder = '../'+folder+'Spectrum_' + str(dt.datetime.now())[:-16] +'/'
else:
	specfolder = '../'+folder+'Spectrum_' + speci_dat+'/'


fig = plt.figure()
max_ener = 0
for spectrum in os.listdir(specfolder):
	# Avoid special files .DS or pdf
	if spectrum[0] == '.' or spectrum[-1] == 'f':
		continue

	ele = float(re.findall('\d+.\d+',spectrum)[0])

	cm_energy_spec = np.load(specfolder+spectrum)

	if max(cm_energy_spec)>max_ener:
		max_ener = max(cm_energy_spec)
	plt.plot(ele*np.ones(len(cm_energy_spec)), cm_energy_spec,'o',color='r')

plt.plot([elemin,elemax],[m0,m0],'-',color='0.5')

if np.dot(d,d)!= 0:
	plt.plot([elemin,elemax],[m1+m2,m1+m2],'-',color='0.5')

def free_part_en(Llat,ntrip):
	conv2 = ((2*np.pi*hbar)/Llat)**2
	n2 = d - ntrip
	Elab = np.sqrt(m1**2+conv2*np.dot(ntrip,ntrip))+np.sqrt(m2**2+conv2*np.dot(n2,n2))
	Ecm = np.sqrt(Elab**2-conv2*np.dot(d,d))
	return Ecm

Ecmsplot = []
for nx in range(4):
	for ny in range(4):
		for nz in range(4):
			ntri = np.array([nx,ny,nz])
			proof =free_part_en(elemin,ntri)
			if proof > max_ener+50:
				continue
			cont = False
			for ecms in Ecmsplot:
				# No two energies closer than 1 MeV (cuz they will be the same)
				if np.abs(proof-ecms) < 1:
					cont = True
			if cont: continue

			Ecmsplot.append(proof)

			Lspace = np.linspace(elemin, elemax, num=100)
			plt.plot(Lspace,free_part_en(Lspace,ntri),'--',color='k')




plt.title(r'Spectrum with BW ($m_0$ = '+str(m0)+' MeV, $g_0$ = '+ str(g0) +')')

plt.xlabel(r'$L/$fm')

plt.ylabel(r'$E/$MeV')

figname = specfolder + 'Spectrum BW (m_0 = '+str(m0)+' MeV, g_0 = '+ str(g0) +')--' + str(dt.datetime.now())[:-16]+'.pdf'

fig.savefig(figname)

plt.show()
