#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def getfRe(fname):
	f = np.loadtxt(fname+'.csv', delimiter=',', skiprows=1)

	p = f[:,0] 
	u = f[:,1:4] 
	assert(np.shape(u)[1] == 3)

	r = f[:,9:]
	assert(np.shape(r)[1] == 3)

	# constants
	ubar = 1 # average velocity 
	D = 50e-3
	d = .2 
	A = D**2
	P = 4*D 
	Dh = 4*A/P 
	nu = (1.004e-3 + 9.8313e-4)/2 # nu = mu/rho 

	# reynolds 
	Re = ubar*Dh/nu 
	print('Re =', Re)

	v = np.sqrt(u[:,0]**2 + u[:,1]**2 + u[:,2]**2) # velocity magnitude 

	# dim v 
	lbound = .2
	ubound = .4

	z1 = np.argmin(np.fabs(r[:,0] - lbound))
	z2 = np.argmin(np.fabs(r[:,0] - ubound))

	dp = p[z2] - p[z1] # dp/rho 
	dz = r[z2,0] - r[z1,0]

	fRe = -dp/dz * Dh/(.5*ubar**2) * ubar*Dh / nu 
	print('fRe =', fRe)

	fanning = 8*-dp*(A/P)**2 * ubar/nu/dz 
	print('fanning =', fanning)	

	return fRe

getfRe('axial')

# f = np.zeros(3)
# dx = np.zeros(3)

# f[0] = getfRe('axial160')
# f[1] = getfRe('axial80')
# f[2] = getfRe('axial40')

# r = 2

# pobs = np.log(np.fabs((f[2] - f[1])/(f[1] - f[0])))/np.log(r)

# print('pobs =', pobs)

# ptrue = 2
# Fs = 1.25 # factor of safety 
# if (np.fabs(ptrue - pobs)/ptrue > .1): # update if convergence off 
# 	Fs = 3 

# # compute GCI
# gci = Fs/(r**pobs - 1)*np.fabs(f[1] - f[0])

# print('fRe =', f[0], '+-', gci)