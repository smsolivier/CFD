#!/usr/bin/env python3 

import numpy as np 
import matplotlib.pyplot as plt 

import os 
import sys 

from scipy.interpolate import interp1d 
from scipy.integrate import quadrature 
import scipy.integrate as integrate 
from scipy.stats import norm

def calcGCI(f1, f2, f3, r12, r23, tol=1e-6):
	''' computes observed order of convergence and sigma from 
		non uniform refinement 
	''' 
	converged = 0 
	pold = 2 
	alpha = np.fabs((f3 - f2)/(f2 - f1)) 
	while (not(converged)):
		p = np.log(alpha* (r12**pold - 1)/(r23**pold - 1))/np.log(r12)

		if (np.fabs(p - pold)/p < tol):
			converged = 1 

		pold = p 

	Fs = 3 
	if (np.fabs(p - 2)/2 < .1):
		Fs = 1.25 

	sigma = Fs/(r12**p - 1) * np.fabs(f1 - f2) 

	return p, sigma 

def gaussian(mu, sigma):
	''' return gaussian with average mu and variance sigma^2 ''' 
	f = lambda x: 1/(np.sqrt(2*sigma**2*np.pi))*np.exp(
		-1/2*((x-mu)/sigma)**2)

	# f = lambda x: 1/(sigma*np.sqrt(2*np.pi))*np.exp(-1/2*((x - mu)/sigma)**2)

	return f 

def calcOmega(u_sim, sigma_sim, u_ex, sigma_ex):
	''' return overlapping fidelity Omega ''' 
	a = 1/(2*sigma_sim*np.sqrt(np.pi))
	b = lambda x: ((x - u_ex)/sigma_ex)**2 
	c = lambda x: ((x - u_sim)/sigma_sim)**2 

	omega = lambda x: a * np.exp(-1/4*(b(x) + c(x)))

	# x = np.linspace(u_sim-2, u_sim+2, 100)
	# plt.plot(x, omega(x))
	# plt.title(str(u_sim/u_ex))
	# plt.show()

	x = np.linspace(u_sim - .5, u_sim + .5, 100)
	f_ex = gaussian(u_ex, sigma_ex)
	f_sim = gaussian(u_sim, sigma_sim)

	# plt.plot(x, f_ex(x), label='ex')
	# plt.plot(x, f_sim(x), label='sim')
	# plt.plot(x, omega(x), label='omega')
	# plt.legend(loc='best')
	# plt.show()

	# print(u_sim/u_ex, sigma_sim/sigma_ex)
	# Omega, err = quadrature(omega, u_sim-1, u_sim+1)

	# integrate omega(x) dx 
	Omega = integrate.quad(omega, -np.inf, np.inf)

	# print('Omega =', Omega, 'Error =', err)

	return Omega[0] 

def readExperiment(case, dist):
	''' read in experimental data for location dist ''' 
	if (case == 0):
		pre = 'N320/'
		post = '_PIV_dr1_u06.dat'
	elif (case == 1):
		pre = 'N337/'
		post = '_PIV_dr0_u10.dat'
	elif (case == 2):
		pre = 'N339/'
		post = '_PIV_dr0_u06.dat'

	pre = 'data/' + pre + 'X'

	# y, Vx, U_x95, k, kLow, kUp
	df = np.loadtxt(pre+dist+post, skiprows=5, usecols=(1,2,3,15,16,17)) 

	y_ex = df[:,0] # y grid points 
	u_ex = -df[:,1] # Ux profile 
	sigma_ex = df[:,2] # 95% uncertainty 
	k = df[:,3]
	kLow = df[:,4]
	kUp = df[:,5] 

	ksigma = (kLow - kUp)/2 

	return y_ex, u_ex, sigma_ex, k, ksigma 

def readGCI(dr, grid, u_ex, sigma_ex, k_ex, ksigma_ex):
	''' read in data from runGCI 
		interpolates onto grid 
		calculates sigma from GCI 
		computes Omega 
		returns first simulation, uncertainty and omega 
			at each point in grid 
	''' 
	N = np.zeros(3) # store number of volumes for each run 

	U = np.zeros((3,len(grid))) # store velocities on grid 
	K = np.zeros((3, len(grid))) # store k for each run on each grid point 
	C = np.zeros((3, len(grid))) # store C for each run on each grid point 
	centerU = np.zeros(3) # store centerline velocity for each run 
	centerk = np.zeros(3) # store centerline k for each run 
	for i in range(3): # loop through number of runs 
		# get number of volumes 
		fname = dr + '/run'+str(i)+'.txt'
		f = open(fname, 'r')
		N[i] = float(f.readline().strip())

		# read in data 
		y, u, k, c = np.loadtxt(fname, skiprows=1, unpack=True)

		# interpolate onto grid 
		# U[i,:] = np.interp(grid, y, u) # interpolate onto grid 
		U[i,:] = interp1d(y, u, kind='cubic')(grid) 
		kfunc = interp1d(y, k, kind='cubic') # interpolated function, used for integration 
		K[i,:] = kfunc(grid) # interpolate k onto grid 
		C[i,:] = interp1d(y, c, kind='cubic')(grid)

		centerU[i] = U[i,int(len(grid)/2)] # centerline Ux 
		centerk[i] = K[i,int(len(grid)/2)] # centerline k 

	# calculate refinement factor 
	r12 = (N[0]/N[1])**(1/2)
	r23 = (N[1]/N[2])**(1/2) 

	pU, sigmaU = calcGCI(centerU[0], centerU[1], centerU[2], r12, r23) 

	pk, sigmak = calcGCI(centerk[0], centerk[1], centerk[2], r12, r23) 

	Omega = calcOmega(centerU[0], sigmaU, u_ex[int(len(grid)/2)], sigma_ex[int(len(grid)/2)])

	return U[0,:], 2*sigmaU*np.ones(len(grid)), Omega*np.ones(len(grid)) 

def globalMerit(Omega, y_ex, u_ex, u, alpha, beta): 
	''' use omega and differnence in derivatives to calculate global merit 
		Small M is better 
		uses backward difference to evaluate derivatives 
			lose first data point 
	''' 

	# create comparison grid 
	# grid = np.linspace(y_ex[0], y_ex[-1], 100)
	grid = y_ex 

	# create experimental data 
	# interpolate exact answer onto grid 
	# ex = interp1d(y_ex, u_ex)(grid)
	ex = u_ex 

	# create simulation data 
	# interpolate simulation onto grid 
	# sim = interp1d(y_ex, u)(grid) 
	sim = u 

	# store derivatives 
	Ndata = len(grid) - 1 # number of evaluated grid points 
	dex = np.zeros(Ndata) # derivative of experimental 
	dsim = np.zeros(Ndata) # derivative of simulation 

	# compute derivatives 
	for i in range(1, len(grid)):
		# backward difference derivatives
		dex[i-1] = (ex[i] - ex[i-1])/(grid[i] - grid[i-1]) 
		dsim[i-1] = (sim[i] - sim[i-1])/(grid[i] - grid[i-1]) 

	# plt.plot(grid[1:], dex)
	# plt.plot(grid[1:], dsim)
	# plt.show()

	E = np.fabs(1 - dsim/dex) # error in derivatives 

	# compute M 
	M = 0
	for i in range(1, len(grid)):
		M += alpha*(1 - Omega[i]) + beta*E[i-1]
	M /= Ndata 

	return M 

def handle(case, expDir, ofDir, alpha, beta):
	''' combines 
			readExperiment 
			readGCI 
			globalMerit
	''' 

	# get experimental data 
	y_ex, u_ex, sigma_ex, k_ex, ksigma_ex = readExperiment(case, expDir)

	# get GCI data 
	u, sigma, Omega = readGCI(ofDir, y_ex/1000, u_ex, sigma_ex, k_ex, ksigma_ex) 

	# compute M 
	M = globalMerit(Omega, y_ex, u_ex, u, alpha, beta)

	return y_ex, u_ex, sigma_ex, u, sigma, Omega, M 

# undisclosed weighting factors for M formula 
alpha = 1 
beta = 1 

case = 2 
dist = ['050', '150', '250', '350', '450'] 
gciDir = 'gciData/'
subDirs = sorted(os.listdir(gciDir)) # get list of subdirectories in gciDir 

for i in range(len(dist)):
	y_ex, u_ex, sigma_ex, u, sigma, Omega, M = handle(
		case = case, 
		expDir = dist[i], 
		ofDir = gciDir + subDirs[i], 
		alpha = alpha, 
		beta = beta)

	print(dist[i], 'M =', M)

	# print(Omega)

	plt.subplot(np.ceil(len(dist)/2), 2, i+1) 
	plt.errorbar(y_ex, u, yerr=sigma)
	plt.plot(y_ex, u_ex)
	plt.title('M = ' + str(M))

	# print('\n', i)
	# for i in range(len(y_ex)):
	# 	print(u[i], sigma[i], u_ex[i], sigma_ex[i])

plt.show()

# print(calcOmega(1, 1e-4, 1, .1))