#!/usr/bin/env python3

import Timer 
import os
import sys 

import numpy as np 

from scipy.interpolate import interp1d
import scipy.integrate as integrate 

### 

# Read in the data in a useable format
def readin(filename, skip=0):
	"""Reads in the output file"""
	f = open(filename, 'r')
	with open(filename) as f:
	    lines = f.readlines()
	f.close()
	return lines[skip:]


def maxVals(data, LOUD=False):
	"""Determines max velocity and TKE with uncertainty of given inlet conditions.
	Return:
		[max velocity], [velocity uncertainty], [max TKE], [TKE uncertainty]
	"""
	data = ''.join(data)
	
	conversion = 1000

	# Reinitializes the data file name
	if (data == '0.6' or data == '.6'):
		data = 'Vel_Inlet_X-50_u06ms.dat'
	elif (data == '1.0' or data == '1'):
		data = 'Vel_Inlet_X-50_u10ms.dat'
	else:
		error('Invalid inlet conditions, either \'0.6\' or \'1.0\'.')

	# Determine the directories
	projectPATH = os.path.dirname(os.path.realpath(__file__))

	#####################################################################
	# Load in the data for interpolation
	#####################################################################
	# print("Reading inlet data from '{}'... ".format(projectPATH+'/data/'+data), end='\n')
	lines = readin(projectPATH+'/data/'+data, skip=5)
	cols = np.array(["x", "y", "z", "Vx", "U_x95", "Vy", "U_y95", "Vz", "U_z95",
		"RMSVx", "RMSVxUpper", "RMSVxLower", "RMSVy", "RMSVyUpper", "RMSVyLower",
		"RMSVz", "RMSVzUpper", "RMSVzLower", "k", "kLOWER", "kUPPER"])
	# Remove the black lines at the end of the file
	for i in range(len(lines))[::-1]:
		if (lines[i].strip() == ''):
			del lines[-1]
		else:
			break
	nRows = len(lines)
	nCols = len(cols)

	# Create a matrix of all the values from the file
	matdata = np.zeros([nRows, nCols])
	for i in range(len(lines)):
		matdata[i,:] = np.fromstring(lines[i], sep='\t')
		# print(i, matdata[i,:]) # Data verification purposes

	# for i in matdata:
	# 	print(i)

	# print(matdata[30])

	##########
	# 'kdata' contians the corresponding turbulent kinetic energy data
	# 'ukdata' holds the lower (FIRST COLUMN) and upper (SECOND COLUMN) bounds of k
	# 'veldata' contains the corresponding velocity component values
	# 'uveldata' holds the standard deviation of velocity in the x, y, and z
	##########
	veldata = np.zeros([nRows, 3])
	uveldata = np.zeros([nRows, 3])
	kdata = np.zeros(nRows)
	ukdata = np.zeros([nRows, 2])
	for i in range(nRows):
		kdata[i] = np.array([matdata[i,18]])
		ukdata[i,:] = np.array([matdata[i,19], matdata[i,20]])
		veldata[i,:] = np.array([matdata[i,3], matdata[i,5], matdata[i,7]])
		# print(np.array([matdata[i,3], matdata[i,5], matdata[i,7]]))
		uveldata[i,:] = np.array([matdata[i,4], matdata[i,6], matdata[i,8]])/2
	#####################################################################

	vmax = 0
	vindex = 0
	for i, row in enumerate(veldata):
		if (row[0] > vmax):
			vmax = row[0]
			vindex = i

	kmax = 0
	kindex = 0
	for i, val in enumerate(kdata):
		if (val > kmax):
			kmax = val
			kindex = i

	if (LOUD == True):
		print("Velocity:", veldata[vindex,0], "\ Uncertainty:", uveldata[vindex,0])
		print("TKE:", kdata[kindex], "\ Uncertainty:", kdata[kindex])

	return veldata[vindex,0], uveldata[vindex,0], kdata[kindex], kdata[kindex]

# sensitivity on velocity, k, mu, rho, D

def get(dr):

	y, k = np.loadtxt(dr+'/k', unpack=True) # read in k from dr directory
	y, U = np.loadtxt(dr+'/U', unpack=True)
	y, C = np.loadtxt(dr+'/C', unpack=True)

	kfunc = interp1d(y, k, kind='cubic') # make interpolated k function
	Ufunc = interp1d(y, U, kind='cubic')
	Cfunc = interp1d(y, C, kind='cubic')

	# integrate k 
	kint = integrate.quad(kfunc, y[0], y[-1])[0]
	Uint = integrate.quad(Ufunc, y[0], y[-1])[0]
	Cint = integrate.quad(Cfunc, y[0], y[-1])[0] 

	return kint, Uint, Cint


readDir = 'output/450' # place to read k values from 

# number of volumes 
Nx1 = 10 
Nx2 = 30 
Ny = 20 
Nz = 1 

# parallel decomposition
npx = 4 
npy = 1
npz = 1

# case to run 
case = 0 

if (case == 0):
	file = 'N339'
elif (case == 1):
	file = 'N337'
elif (case == 2):
	file = 'N320'
else:
	file = 'N318'

if len(sys.argv) == 2:
	case = int(sys.argv[1])

# set case dependent values 
inletVelocity = '' # store inlet velocity string for generateInlet  
if (case == 0): # N339, .6 m/s, same density 
	inletVelocity = '.6' 

	# set density 
	rhoTop = 998 
	rhoBottom = 998 

	# set mu 
	muTop = 1.002e-3 
	muBottom = 1.002e-3 

elif (case == 1): # N337, 1 m/s, same density  
	inletVelocity = '1' 

	# set density 
	rhoTop = 998
	rhoBottom = 998 

	# set mu 
	muTop = 1.002e-3 
	muBottom = 1.002e-3 

elif (case == 2): # N320, .6 m/s, different density  
	inletVelocity = '.6'

	# set density 
	rhoTop = 998 
	rhoBottom = 1008 

	# set mu 
	muTop = 1.002e-3 
	muBottom = .991e-3 

elif (case == 3): # N318, 1 m/s, different density, blind case 
	inletVelocity = '1'

	# set density 
	rhoTop = 998 
	rhoBottom = 1008 

	# set mu 
	muTop = 1.002e-3 
	muBottom = .991e-3 

else: # exit if other case given 
	print('no case found')
	sys.exit()

# set mass diffusivity 
Dab = 1e-6 

# number of simulations per variable 
nruns = 5 

# multiplicative factor for each variable 
scale = np.logspace(-3, -.5, nruns)

kint = np.zeros(nruns)
Uint = np.zeros(nruns)
Cint = np.zeros(nruns)

# velocity, k, mu, rho, D
coefvals = np.zeros([3,7])

# get k with no perturbation 
# run openfoam 
runstr = './run \
	-N {} {} {} {} \
	-np {} {} {} \
	-case {} \
	-quiet \
	>senslog'.format(
		Nx1, Nx2, Ny, Nz, # number of volumes 
		npx, npy, npz, # parallel decomposition 
		case, # case to run 
		)

x = os.system(runstr)
if (x != 0): # exit if problems 
	print(runstr)
	sys.exit()
k_0, U_0, C_0 = get(readDir)

vmax, uv, kmax, uk = maxVals(inletVelocity)
namespace = "Scale,k,U,C"

print("Ux Sensitivity")
coef = np.zeros([3,nruns])
for i in range(nruns):
	print("Scale {}:".format(scale[i]))
	# run openfoam 
	runstr = './run \
		-N {} {} {} {} \
		-np {} {} {} \
		-case {} \
		-Ufactor {} \
		-quiet \
		>senslog'.format(
			Nx1, Nx2, Ny, Nz, # number of volumes 
			npx, npy, npz, # parallel decomposition 
			case, # case to run 
			scale[i] # scale velocity 
			)

	x = os.system(runstr)
	if (x != 0): # exit if problems 
		print(runstr)
		sys.exit()

	k, U, C = get(readDir)
	coef[0,i] = (k - k_0)/(vmax*(scale[i]-1))
	coef[1,i] = (U - U_0)/(vmax*(scale[i]-1))
	coef[2,i] = (C - C_0)/(vmax*(scale[i]-1))
	print("\tk:", coef[0,i])
	print("\tU:", coef[1,i])
	print("\tC:", coef[2,i])

for i, row in enumerate(coef):
	coefvals[i,:] = np.mean(row)
	if (np.std(row) > 0.5*np.mean(row)):
		if (i == 0):
			print("Error: large standard deviation in k.")
		if (i == 1):
			print("Error: large standard deviation in U.")
		if (i == 2):
			print("Error: large standard deviation in C.")
temp = np.concatenate([scale.reshape([1,nruns]),coef], axis=0)
temp = temp.T
np.savetxt('InputUncertainty/Runs/{}-U.csv'.format(file), temp, header=namespace, delimiter=',')

print("k Sensativity")
coef = np.zeros([3,nruns])
for i in range(nruns):
	print("Scale {}:".format(scale[i]))
	# run openfoam 
	runstr = './run \
		-N {} {} {} {} \
		-np {} {} {} \
		-case {} \
		-kfactor {} \
		-quiet \
		>senslog'.format(
			Nx1, Nx2, Ny, Nz, # number of volumes 
			npx, npy, npz, # parallel decomposition 
			case, # case to run 
			scale[i] # scale velocity 
			)

	x = os.system(runstr)
	if (x != 0): # exit if problems 
		print(runstr)
		sys.exit()

	k, U, C = get(readDir)
	coef[0,i] = (k - k_0)/(kmax*(scale[i]-1))
	coef[1,i] = (U - U_0)/(kmax*(scale[i]-1))
	coef[2,i] = (C - C_0)/(kmax*(scale[i]-1))
	print("\tk:", coef[0,i])
	print("\tU:", coef[1,i])
	print("\tC:", coef[2,i])

for i, row in enumerate(coef):
	coefvals[i,:] = np.mean(row)
	if (np.std(row) > 0.5*np.mean(row)):
		if (i == 0):
			print("Error: large standard deviation in k.")
		if (i == 1):
			print("Error: large standard deviation in U.")
		if (i == 2):
			print("Error: large standard deviation in C.")
temp = np.concatenate([scale.reshape([1,nruns]),coef], axis=0)
temp = temp.T
np.savetxt('InputUncertainty/Runs/{}-k.csv'.format(file), temp, header=namespace, delimiter=',')

coef = np.zeros([3,nruns])

def run(dr, Nx1, Nx2, Ny, Nz, npx, npy, npz, case, props):
	muTop = props[0]
	rhoTop = props[1]
	muBottom = props[2]
	rhoBottom = props[3]
	Dab = props[4]

	runstr = './run \
		-N {} {} {} {} \
		-np {} {} {} \
		-case {} \
		-muTop {} \
		-rhoTop {} \
		-muBottom {} \
		-rhoBottom {} \
		-Dab {} \
		-quiet \
		>senslog'.format(
			Nx1, Nx2, Ny, Nz, # number of volumes 
			npx, npy, npz, # parallel decomposition 
			case, # case to run 
			muTop, 
			rhoTop, 
			muBottom, 
			rhoBottom,
			Dab
			)

	x = os.system(runstr)
	if (x != 0): # exit if problems 
		print(runstr)
		sys.exit()

	return get(dr)

# properties; muTop, rhoTop, muBottom, rhoBottom, Dab 
propName = np.array(['muTop', 'rhoTop', 'muBottom', 'rhoBottom', 'Dab'])
props = np.array([muTop, rhoTop, muBottom, rhoBottom, Dab])
for i in range(len(props)): # loop through properties
	print("{} Sensativity".format(propName[i]))
	for j in range(nruns): # loop through perturbations
		print("Scale {}:".format(scale[i]))
		newprop = np.copy(props)
		newprop[i] += scale[j]*props[i] # perturb each variable 

		k, U, C = get(readDir)
		coef[0,i] = (k - k_0)/(newprop[i] - props[i])
		coef[1,i] = (U - U_0)/(newprop[i] - props[i])
		coef[2,i] = (C - C_0)/(newprop[i] - props[i])
		print("\tk:", coef[0,i])
		print("\tU:", coef[1,i])
		print("\tC:", coef[2,i])

	for n, row in enumerate(coef):
		coefvals[n,:] = np.mean(row)
		if (np.std(row) > 0.5*np.mean(row)):
			if (n == 0):
				print("Error: large standard deviation in k.")
			if (n == 1):
				print("Error: large standard deviation in U.")
			if (n == 2):
				print("Error: large standard deviation in C.")

	temp = np.concatenate([scale.reshape([1,nruns]),coef], axis=0)
	temp = temp.T
	np.savetxt('InputUncertainty/Runs/{}-{}'.format(file,propName[i]), temp, header=namespace, delimiter=',')
	coef = np.zeros([3,nruns])

inputs = np.array(['Ux', 'k', 'mu_top', 'mu_bottom', 'rho_top', 'rho_bottom', 'Dab'])
# uncertainty[-1] is the standard deviation of Dab, assumed to be 5e-7 with a value
# of 1e-6 (round off error)
uncertainty = np.array([uv, uk, 0.0005, 0.0005, 0.5, 0.5, 5e-7])
inputu = 0
fmt='%.4e'

sensativities = np.array(['k', 'U', 'C'])

for n, sen in enumerate(sensativities):
	f = open('InputUncertainty/{}-{}.csv'.format(file,sen), 'w')
	f.write('Input,Coefficient,Uncertainty\n')
	for i, j, k in zip(inputs, coefvals[n,:], uncertainty):
		f.write('{},{},{}\n'.format(i,fmt%j,fmt%k))
		inputu += (j*k)**2
	f.write('Overall Uncertainty,{}\n'.format(np.sqrt(inputu)))
	f.close()
