#!/usr/bin/env python3

import Timer 
import os
import sys 

import numpy as np 

from scipy.interpolate import interp1d
import scipy.integrate as integrate 

# sensitivity on velocity, k, mu, rho, D

def get(dr, Nx1, Nx2, Ny, Nz, npx, npy, npz, case, ufactor):

	# run openfoam 
	runstr = './run \
		-N {} {} {} {} \
		-np {} {} {} \
		-case {} \
		-Ufactor {} \
		>senslog'.format(
			Nx1, Nx2, Ny, Nz, # number of volumes 
			npx, npy, npz, # parallel decomposition 
			case, # case to run 
			ufactor # scale velocity 
			)

	x = os.system(runstr)
	if (x != 0): # exit if problems 
		print(runstr)
		sys.exit()

	y, k = np.loadtxt(dr+'/k', unpack=True) # read in k from dr directory 

	kfunc = interp1d(y, k, kind='cubic') # make interpolated k function 

	# integrate k 
	kint = integrate.quad(kfunc, y[0], y[-1])[0] 

	return kint 

readDir = 'output/350' # place to read k values from 

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

# number of simulations per variable 
nruns = 5 

# multiplicative factor for each variable 
scale = np.logspace(-3, -.5, nruns)

kint = np.zeros(nruns)

# get k with no perturbation 
k_0 = get(readDir, Nx1, Nx2, Ny, Nz, npx, npy, npz, case, 1)

# store coefficient from each run 
coef = np.zeros(nruns)
for i in range(nruns):
	k = get(readDir, Nx1, Nx2, Ny, Nz, npx, npy, npz, case, scale[i])

	coef[i] = (k - k_0)/scale[i]

	print(coef[i])

print(np.mean(coef), np.std(coef))

