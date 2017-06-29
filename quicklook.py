import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

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
		if line[0:2] == 'sd':
			sd = True
			speci_dat = str(line[3:-1])
		if line[0:2] == 'g0':
			g0 = float(line[2:]) # coupling
		if line[0:2] == 'M0':
			m0 = float(line[3:]) # MeV, Resonance mass

d = [dx,dy,dz]

# -----------------
# Define constants
I = complex(0.0,1.0)
hbar = 197.3269718 # MeV fm when speed of ligth is 1

# -----------------
# Define variables of the Z function
# d= [1.0,1.0,0]
# L_lattice_size = 10.0 # fermis
L = L_lattice_size/hbar # MeV^-1
lab_moment = (2*np.pi/L)*np.array(d)
lab_moment2 = np.dot(lab_moment,lab_moment)
m1 = 146.0 # MeV
m2 = 146.0 # MeV
# m0 = 770.0 # MeV, Resonance mass
# g0 = 3.0 # coupling

# save_folder = '../'+folder+'Plots/'

filename = 'Z_' + str(l) +'_' + str(m)+ '_d(' + str(d) +')'

# For plotting a specific data set use the following format:
# folder = 'Z_values_'+'2017-mm-dd'+'/'

if sd == False:
	zfolder = '../'+folder+'Z_values_'+str(dt.datetime.now())[:-16]+'/'
else:
	zfolder = '../'+folder+'Z_values_'+speci_dat+'/'

zdlm = np.load(zfolder+filename+'.npy')

filename = 'x2-'+filename 

x2 = np.load(zfolder+filename+'.npy')

q2 = x2*((2*np.pi)/L)**2

def cdlm(l,x2,zdlm):
	# Kinematics
	cm_energy = np.sqrt(q2 + m1**2) + np.sqrt(q2 + m2**2)
	lab_energy = np.sqrt(cm_energy**2+lab_moment2)
	gamma = lab_energy/cm_energy

	firs_fact = np.sqrt(4*np.pi)/(gamma * L**3)*(2*np.pi/L)**(l-2)

	value = firs_fact*zdlm

	return value

def pcotd(q2):
	# Kinematics
	cm_energy = np.sqrt(q2 + m1**2) + np.sqrt(q2 + m2**2)

	val = (m0**2-cm_energy**2)*6*np.pi*cm_energy/(g0*m0)**2

	return val

def diff(q2):
	x2 = q2*(L/(2*np.pi))**2
	return np.abs(np.real(4*np.pi*cdlm(0,x2,zdlm))-pcotd(q2))

diff_arr = diff(q2)

# Guess Luscher solutions
q2mask = np.ma.masked_array(q2,np.isnan(diff_arr))
diff_arr_mask = np.ma.masked_array(diff_arr,np.isnan(diff_arr))

q2_2gues = q2mask.compressed()
diff_arr2gues = diff_arr_mask.compressed()

guesses = np.array([])
g1 = diff_arr2gues[1]

for ii in range(len(diff_arr2gues)-2):
	if g1 < diff_arr2gues[ii-1] and  g1 < diff_arr2gues[ii+1]:
		guesses = np.append(guesses,q2_2gues[ii])

	g1 = diff_arr2gues[ii+1]

num_gues =len(guesses)

print 'Guesses found: ' + str(num_gues)

# Guess annotation
vert_pos = .88

guestring = 'Guess at: \n $p^{2} \in \{$'
for ii in xrange(num_gues):
	if ii+1 == num_gues:
		guestring += '%.2f' % guesses[ii] + '$\}$'
	elif ii < 5 and (ii+1)%5 == 0 :
		vert_pos -= 0.03
		guestring += '%.2f' % guesses[ii] +',\n'
	elif (ii+1)%6  == 0 and ii > 5:
		vert_pos -= 0.03
		guestring += '%.2f' % guesses[ii] +',\n'
	else:
		guestring += '%.2f' % guesses[ii] +', '	

# ---------
# Matplotlib param
plt.rc('text', usetex=True)
plt.rc('font', size=12)
plt.rc('font', family='serif')
plt.rc('axes.formatter', useoffset = False)

# Plot
fig = plt.figure()
plt.plot(q2,np.real(4*np.pi*cdlm(0,x2,zdlm)),'-',q2,pcotd(q2),'-',q2,diff_arr,'-')

# tidy up the figure
plt.grid(True)
plt.xlabel(r'$p^2$')
plt.legend(['$4\pi c_{00}^{[000]}$', '$p\cot\delta_{S}$', 'difference'], loc='lower right')
plt.annotate(guestring, xy=(0.05, vert_pos), xycoords="axes fraction",
			bbox=dict(boxstyle="round",fc="1"))

#plt.show()

save_folder = '../'+folder+'Guess_E_' + str(dt.datetime.now())[:-16] +'/'

if not os.path.exists(save_folder):
	os.makedirs(save_folder)

filename = 'L_'+ str(L_lattice_size) +'_'+ str(l) +'_' + str(m)+ '_d(' + str(d) +')'

figname = filename +'.pdf'

fig.savefig(save_folder + figname)

np.save(save_folder+filename,guesses)
