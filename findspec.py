import os
import sys
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import scipy.optimize as opt
import datetime as dt
import numpy as np

# -----------------
# Config file
config_file = str(sys.argv[1])
sd = False

with open(config_file, 'r') as f:
	for line in f:
		if line[0:7] == 'folder ':
			folder = line[7:-1]
		if line[0:2] == 'l ':
			l = int(line[2:])
		if line[0:2] == 'm ':
			m = int(line[2:])
		if line[0:3] == 'dx ':
			dx = int(line[2:])
		if line[0:3] == 'dy ':
			dy = int(line[2:])
		if line[0:3] == 'dz ':
			dz = int(line[2:])
		if line[0:2] == 'L ':
			L_lattice_size = float(line[2:]) # fermis
		if line[0:6] == 'x2max ':
			x2max = float(line[6:])
		if line[0:6] == 'x2min ':
			x2min = float(line[6:])
		if line[0:7] == 'points ':
			points = int(line[7:])
		if line[0:3] == 'sd ':
			sd = True
			speci_dat = str(line[3:-1])
		if line[0:3] == 'g0 ':
			g0 = float(line[3:])
		if line[0:3] == 'M0 ':
			m0 = float(line[3:])
			
d = [dx,dy,dz]

# -----------------
# Define constants
I = complex(0.0,1.0)
hbar = 197.3269718 # MeV fm when speed of ligth is 1

# Define variables of the C function
# S-wave
# d= [0.0,0.0,0]
# L_lattice_size = 10.0 # fermis
L = L_lattice_size/hbar # MeV^-1
lab_moment = (2*np.pi/L)*np.array(d)
lab_moment2 = np.dot(lab_moment,lab_moment)
m1 = 146.0 # MeV
m2 = 146.0 # MeV
# m0 = 770.0 # MeV, Resonance mass
# g0 = 3.0 # coupling


# Get guesses
if sd == False:
	guessfolder = '../'+folder+'Guess_E_'+str(dt.datetime.now())[:-16]+'/'
else:
	guessfolder = '../'+folder+'Guess_E_'+speci_dat+'/'

filename = 'L_'+ str(L_lattice_size) +'_'+ str(l) +'_' + str(m)+ '_d(' + str(d) +')'

guesses = np.load(guessfolder+filename+'.npy')

x2spacing = (x2max-x2min)/(points-1.0)

def cfdlm(l,q2):
	# Kinematics
	cm_energy = np.sqrt(q2 + m1**2) + np.sqrt(q2 + m2**2)
	lab_energy = np.sqrt(cm_energy**2+lab_moment2)
	gamma = lab_energy/cm_energy

	firs_fact = np.sqrt(4*np.pi)/(gamma * L**3)*(2*np.pi/L)**(l-2)

	x2 = q2*(L/(2*np.pi))**2

	value = firs_fact*fdlm(x2)

	return value

def pcotdf(q2):
	# Kinematics
	cm_energy = np.sqrt(q2 + m1**2) + np.sqrt(q2 + m2**2)
	val = (m0**2-cm_energy**2)*6*np.pi*cm_energy/(g0*m0)**2

	# Poles
	q2_arr = np.repeat([q2],num_poles,axis=0)
	pol_arr = np.transpose(np.repeat([poles],len(q2),axis=0))

	factor = q2_arr*(L/(2*np.pi))**2-pol_arr
	gggg= np.prod(factor,axis=0)

	return gggg*val

cm_energy_spec = []

save_folder = '../'+folder+'Solution_E_' + str(dt.datetime.now())[:-16] +'/'

if not os.path.exists(save_folder):
	os.makedirs(save_folder)

# Get the F values and solve for E
for guess in guesses:

	if sd == False:
		ffolder = '../'+folder+'F_values_g_'+str(dt.datetime.now())[:-16]+'/'
	else:
		ffolder = '../'+folder+'F_values_g_'+speci_dat+'/'

	filename = 'F_' + str(l) +'_' + str(m)+ '_d(' + str(d) +')_guess(%.1f)'%guess

	y = np.load(ffolder+filename+'.npy')

	filename = 'x2-'+filename 

	x = np.load(ffolder+filename+'.npy')

	filename = 'poles-'+ filename[3:]

	poles = np.load(ffolder+filename+'.npy')
	num_poles = len(poles)

	xmask = np.ma.masked_array(x,np.isnan(y))
	ymask = np.ma.masked_array(y,np.isnan(y))

	xtoint = xmask.compressed()

	ytoint = ymask.compressed()

	fdlm = interp1d(xtoint, np.real(ytoint), kind='cubic')

	# Plot
	# Matplotlib param
	plt.rc('text', usetex=True)
	plt.rc('font', size=12)
	plt.rc('font', family='serif')
	plt.rc('axes.formatter', useoffset = False)

	# x2plmax = (guess)*(L/(2*np.pi))**2 + 2 * x2spacing
	# x2plmin = (guess)*(L/(2*np.pi))**2 - 2 * x2spacing

	x2plmin = xtoint[0]
	x2plmax = xtoint[-1]
	
	# Ensure plot is inside the interpolation range
	absdif = x2plmax-x2plmin
	x2plmin += .01*absdif
	x2plmax -= .01*absdif

	q2min = x2plmin*((2*np.pi)/L)**2
	q2max = x2plmax*((2*np.pi)/L)**2
	q2 = np.linspace(q2min, q2max, num=1000)

	fig = plt.figure()
	plt.plot(q2, 4*np.pi*cfdlm(0,q2),'-',q2,pcotdf(q2),'-')

	# tidy up the figure
	plt.grid(True)
	plt.xlabel(r'$p^2$')
	cname = 'c_{' + str(l) + str(m)+ '}^{' + str(d) +'}'
	plt.legend(['$4\pi '+ cname + '\prod (q^2-q_0)$', '$p\cot\delta_{S}\prod (p^2-p_0)$'], loc='lower right')
	#plt.show()

	filename = 'L_'+ str(L_lattice_size) +'_'+ str(l) +'_' + str(m)+ '_d(' + str(d) +')_guess(%.1f)'%guess

	figname = filename +'.pdf'

	fig.savefig(save_folder + figname)

	def lus(q2):
		return pcotdf(q2)-4*np.pi*cfdlm(0,q2)

	q2spec = opt.fsolve(lus,guess)
	ener = np.sqrt(q2spec + m1**2) + np.sqrt(q2spec + m2**2)
	cm_energy_spec.append(ener)
	print('Solution found at %.3f MeV' % ener) 


# for nn in range(6):
# 	free_part_en = lambda Llat: 2*np.sqrt((nn*(2*np.pi*hbar)/Llat)**2 + m1**2)
# 	Lspace = np.linspace(3, 12, num=100)
# 	plt.plot(Lspace,free_part_en(Lspace),'--')

# plt.plot(L_lattice_size*np.ones(len(cm_energy_spec)), cm_energy_spec,'o')

# plt.show()

filename = 'Spectrum_L_' + str(L_lattice_size) +'_l'+ str(l) +'_' + str(m)+ '_d(' + str(d) +')'

specfolder = '../'+folder+'Spectrum_' + str(dt.datetime.now())[:-16] +'/'

if not os.path.exists(specfolder):
	os.makedirs(specfolder)

np.save(specfolder+filename,np.array(cm_energy_spec))

