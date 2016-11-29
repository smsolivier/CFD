#!/usr/bin/env python3 

import numpy as np 
import matplotlib.pyplot as plt 

import os 
import sys 

from scipy.interpolate import interp1d 
from scipy.integrate import quadrature 
import scipy.integrate as integrate 
from scipy.stats import norm

def gaussian(mu, sigma):
	''' return gaussian with average mu and variance sigma^2 ''' 
	
	f = lambda x: 1/np.sqrt(2*sigma**2*np.pi) * np.exp(
		-1*(x - mu)**2/(2*sigma**2))


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

	# x = np.linspace(u_sim - .5, u_sim + .5, 100)
	# f_ex = gaussian(u_ex, sigma_ex/2)
	# f_sim = gaussian(u_sim, sigma_sim)

	x = np.linspace(norm.ppf(.01, u_sim, sigma_sim), 
		norm.ppf(.99, u_sim, sigma_sim), 100)
	f_ex = lambda x: norm.pdf(x, u_ex, sigma_ex/2) 
	f_sim = lambda x: norm.pdf(x, u_sim, sigma_sim)

	# print(integrate.quad(f_ex, -np.inf, np.inf)) 
	# print(u_ex, sigma_ex)

	# plt.figure()
	# plt.axvline(u_ex)
	# plt.axvline(u_sim)
	# plt.plot(x, f_ex(x), label='ex')
	# plt.plot(x, f_sim(x), label='sim')
	# plt.plot(x, omega(x), label='omega')
	# plt.legend(loc='best')
	# plt.show()

	# print(u_sim/u_ex, sigma_sim/sigma_ex)
	# Omega, err = quadrature(omega, u_sim-1, u_sim+1)

	# integrate omega(x) dx 
	Omega = integrate.quad(omega, -np.inf, np.inf)

	# print('Omega =', Omega[0], 'Error =', Omega[1])

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

def calcGCI(f, N, tol=1e-6):
	''' computes observed order of convergence and sigma from 
		non uniform refinement 
	'''

	f1 = f[0]
	f2 = f[1]
	f3 = f[2] 

	# check for oscillatory convergence 
	if (f2/f1  - 1> 0 and f2/f3  - 1> 0) or (f2/f1  - 1< 0 and f2/f3  - 1< 0):
		print('oscillatory convergence', f2/f1, f2/f3)

	# check if f increasing with N 
	elif (f1 > f2 > f3):
		print('Increasing trend')

	# asymptotic convergence, good for GCI  
	elif (f1 < f2 < f3):
		print('ok')

	else:
		print('unknown behavior')

	# calculate refinement factor 
	r12 = (N[0]/N[1])**(1/3)
	r23 = (N[1]/N[2])**(1/3) 

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

def readGCI(dr, grid):
	''' read in data from runGCI 
		interpolates onto grid 
		calculates sigma from GCI 
		computes Omega 
		returns first simulation, uncertainty and omega 
			at each point in grid 
	''' 
	N = np.zeros(3) # store number of volumes for each run 

	# Ufunc_ex = interp1d(grid, u_ex, kind='cubic') # exact interpolated velocity 
	# centerU_ex = Ufunc_ex(0) 

	cent = int(len(grid)/2) # center index of grid 

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
		f.close()

		# read in data 
		y, u, k, c = np.loadtxt(fname, skiprows=1, unpack=True)

		# interpolate onto grid 
		U[i,:] = np.interp(grid, y, u) # interpolate onto grid 
		Ufunc = interp1d(y, u, kind='cubic') # interpolated function 
		# U[i,:] = Ufunc(grid) 
		kfunc = interp1d(y, k, kind='cubic') # interpolated function, used for integration 
		K[i,:] = kfunc(grid) # interpolate k onto grid 
		C[i,:] = interp1d(y, c, kind='cubic')(grid)

		centerU[i] = U[i,cent] # centerline Ux 
		centerk[i] = K[i,cent] # centerline k 

	pU, sigmaU = calcGCI(centerU, N) 

	print('pU =', pU, 'U =', centerU[0], sigmaU)

	pk, sigmak = calcGCI(centerk, N) 

	print('pk =', pk, 'k =', centerk[0], sigmak)

	return U[0,:], 2*sigmaU*np.ones(len(grid)), K[0,:], 2*sigmak*np.ones(len(grid)) 

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
	u, sigma, k, sigmak = readGCI(ofDir, y_ex/1000) 

	# compute Omega 
	Omega = np.zeros(len(u))
	for i in range(len(u)):
		Omega[i] = calcOmega(u[i], sigma[i], u_ex[i], sigma_ex[i]) 

	# compute M 
	M = globalMerit(Omega, y_ex, u_ex, u, alpha, beta)

	return y_ex, u_ex, sigma_ex, u, sigma, k_ex, ksigma_ex, k, sigmak, Omega, M 

# undisclosed weighting factors for M formula 
alpha = 1 
beta = 1 

case = 0
dist = ['050', '150', '250', '350', '450'] 
# dist = ['050']
gciDir = 'gciData/'
# gciDir = 'adaData/'
subDirs = sorted(os.listdir(gciDir)) # get list of subdirectories in gciDir 

fig1 = plt.figure()
fig2 = plt.figure()
for i in range(len(dist)):
	y_ex, u_ex, sigma_ex, u, sigma, k_ex, ksigma_ex, k, sigmak, Omega, M = handle(
		case = case, 
		expDir = dist[i], 
		ofDir = gciDir + subDirs[i], 
		alpha = alpha, 
		beta = beta)

	print(dist[i], 'M =', M)

	# print(Omega)

	ax1 = fig1.add_subplot(np.ceil(len(dist)/2), 2, i+1) 
	ax1.errorbar(y_ex, u, yerr=sigma)
	ax1.errorbar(y_ex, u_ex, yerr=sigma_ex)
	ax1.set_title('M = ' + str(M))

	ax2 = fig2.add_subplot(np.ceil(len(dist)/2), 2, i+1)
	ax2.errorbar(y_ex, k, yerr=sigmak)
	# ax2.plot(y_ex, k)
	ax2.errorbar(y_ex, k_ex, yerr=ksigma_ex)

plt.show()

# print(calcOmega(1, 1e-4, 1, .1))