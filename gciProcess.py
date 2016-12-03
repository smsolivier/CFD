#!/usr/bin/env python3 

import numpy as np 
import matplotlib.pyplot as plt 

import os 
import sys 

from scipy.interpolate import interp1d 
from scipy.integrate import quadrature 
import scipy.integrate as integrate 
from scipy.stats import norm

from readFoam import readExperiment # read in experimental values from data/ 

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

def calcGCI(f, N, tol=1e-9):
	''' computes observed order of convergence and sigma from 
		non uniform refinement 
	'''

	f1 = f[0]
	f2 = f[1]
	f3 = f[2] 

	BAD = ''

	# check for oscillatory convergence 
	if (f2/f1  - 1> 0 and f2/f3  - 1> 0) or (f2/f1  - 1< 0 and f2/f3  - 1< 0):
		# print('oscillatory convergence', f2/f1, f2/f3)

		BAD += 'oscillatory ' 
	# asymptotic convergence, good for GCI  
	elif (f1 < f2 < f3 or f1 > f2 > f3):
		BAD += ''

	else:
		# print('unknown behavior')
		BAD += 'unknown '

	# calculate refinement factor 
	r12 = (N[0]/N[1])**(1/3)
	r23 = (N[1]/N[2])**(1/3) 

	converged = 0 # store if converged, break if 1 
	pold = 2 # guess 2 
	alpha = np.fabs((f3 - f2)/(f2 - f1)) # ratio of f's, common for iteration 
	while (not(converged)):
		# update p 
		p = np.log(alpha * (r12**pold - 1)/(r23**pold - 1))/np.log(r12)

		# compare to old p 
		if (np.fabs(p - pold)/p < tol):
			converged = 1 

		# update pold 
		pold = p 

	# set factor of safety 
	Fs = 3 
	if (np.fabs(p - 2)/2 < .1): # if close to expected p = 2
		Fs = 1.25 

	# deviation 
	sigma = Fs/(r12**p - 1) * np.fabs(f1 - f2) 

	if (p < 0): # negative order of convergence 
		# print('p negative')
		BAD += 'negative p '

	if (sigma < 0): # negative uncertainty 
		BAD += 'negative sigma '

	if (p > 10): # p too large 
		# print('p > 20')
		BAD += 'p too large ' 

	return p, sigma, BAD  

def readGCI(dr, grid, cgrid):
	''' read in data from runGCI 
		interpolates onto grid 
		computes metrics for GCI 
		runs GCI on each metric 
			throws out metrics that are bad 
				oscillatory convergence 
				negative p 
				negative sigma 
		uses maximum relative error of surviving metrics 
		applies uniform uncertainty to all grid points 
	''' 
	N = np.zeros(3) # store number of volumes for each run 

	# Ufunc_ex = interp1d(grid, u_ex, kind='cubic') # exact interpolated velocity 
	# centerU_ex = Ufunc_ex(0) 

	cent = int(len(grid)/2) # center index of grid 

	U = np.zeros((3,len(grid))) # store velocities on grid 
	K = np.zeros((3, len(grid))) # store k for each run on each grid point 
	C = np.zeros((3, len(cgrid))) # store C for each run on each grid point 

	# store metrics of interest: centerline Ux, centerline k, integral k, integral Ux
	metricNames = ['U', 'k', 'C', 'kint', 'Uint']
	metrics = np.zeros((len(metricNames), 3)) 
	for i in range(3): # loop through number of runs 
		# get number of volumes 
		fname = dr + '/run'+str(i)+'.txt'
		f = open(fname, 'r')
		N[i] = float(f.readline().strip())
		f.close()

		# read in data 
		y, u, k, c = np.loadtxt(fname, skiprows=1, unpack=True)

		# interpolate onto grid 
		Ufunc = interp1d(y, u, kind='cubic') # interpolated function 
		U[i,:] = Ufunc(grid) # evaluated on grid 
		kfunc = interp1d(y, k, kind='cubic') # interpolated function, used for integration 
		K[i,:] = kfunc(grid) # interpolate k onto grid 
		Cfunc = interp1d(y, c, kind='cubic') # concentration interpolated function 
		C[i,:] = Cfunc(cgrid) # interpolate onto concentration grid 

		# compute GCI metrics 
		metrics[0,i] = Ufunc(0) # interpolated centerline Ux 
		metrics[1,i] = kfunc(0) # interpolated centerline k 
		metrics[2,i] = Cfunc(0) # interpolated centerline concentration 
		metrics[3,i] = integrate.quad(kfunc, grid[0], grid[-1])[0] # integral of k 
		metrics[4,i] = integrate.quad(Ufunc, grid[0], grid[-1])[0] # integral of U 

	# find order of convergence of all metrics, eliminate bad metrics, 
	# find max error of good metrics, use relative error on all return metrics 
	p = np.zeros(len(metrics)) # store order of convergence for each metric 
	sigma = np.zeros(len(metrics)) # store uncertainty for each metric 
	relError = [] # store relative error of good metrics 
	useMetric = [] # store name of good metrics  
	useP = [] # store order of convergence of good metrics 
	for i in range(len(metrics)): # loop through the metrics   
		# get order of convergence 
		p[i], sigma[i], bad = calcGCI(metrics[i,:], N) 

		if (bad == ''): # if a good metric 
			relError.append(sigma[i]/metrics[i,0]) # store relative error of best simulation 
			useMetric.append(metricNames[i]) # store name of good metric 
			useP.append(p[i]) # store convergence of good metric 
		else: # if a bad metric 
			print(metricNames[i] + ': ' + bad) # print reason 

	if (len(relError) == 0): # all metrics bad 
		print('WARNING: no good metrics found')

		error = .1 # set error so that plots ok 

	else: # if at least one good metric 

		relError = np.array(relError) # switch to numpy 
		ind = np.argmax(relError) # location of max error

		error = relError[ind] # maximum relative error 

		print(useMetric[ind], useP[ind], error)

	return ( # parenthesis make able to return on multiple lines 
		U[0,:], U[0,:]*error*np.ones(len(grid)), # return Ux and error 
		K[0,:], K[0,:]*error*np.ones(len(grid)), # return k and error 
		C[0,:], C[0,:]*error*np.ones(len(cgrid)) # return concentration and error 
		)	

def globalMerit(y_ex, u_ex, sigma_ex, u, sigma, alpha, beta): 
	''' use omega and differnence in derivatives to calculate global merit 
		Small M is better 
		uses backward difference to evaluate derivatives 
			lose first data point 
	''' 
	# compute Omega 
	Omega = np.zeros(len(u))
	for i in range(len(u)):
		Omega[i] = calcOmega(u[i], sigma[i], u_ex[i], sigma_ex[i]) 

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

def handle(case, expDir, gciDir, subDir, alpha, beta):
	''' combines 
			readExperiment 
			readGCI 
			globalMerit
	''' 

	ofDir = gciDir + subDir 

	# get experimental data 
	(y_ex, u_ex, sigma_ex, # Ux
		k_ex, ksigma_ex, # k 
		yc, c_ex, sigmac_ex, # concentration 
		x # x location of y profiles 
		) = readExperiment(case, expDir)

	# get GCI data 
	u, sigma, k, sigmak, c, sigmac = readGCI(ofDir, y_ex/1000, yc/1000) 

	# compute M 
	M = globalMerit(y_ex, u_ex, sigma_ex, u, sigma, alpha, beta)

	# compute M for k 
	Mk = globalMerit(y_ex, k_ex, ksigma_ex, k, sigmak, alpha, beta)

	# compute M for C 
	Mc = globalMerit(yc, c_ex, sigmac_ex, c, sigmac, alpha, beta)

	return (
		y_ex, u_ex, sigma_ex, # experemintal Ux 
		u, sigma, # simulated Ux 
		k_ex, ksigma_ex, # experimental k 
		k, sigmak, # simulated k 
		yc, c_ex, sigmac_ex, # experimental concentration 
		c, sigmac, # simulated concentration 
		M, # metric array for Ux 
		Mk, # metric array for k 
		Mc # metric for concentration 
		)  

if __name__ == '__main__':
	# get command line arguments 
	# undisclosed weighting factors for M formula 
	alpha = 1 
	beta = 1 

	# gciDir = 'gciData/'
	gciDir = 'adaData/'

	dist = ['050', '150', '250', '350', '450'] 
	# dist = ['050']
	subDirs = sorted(os.listdir(gciDir)) # get list of subdirectories in gciDir 

	# get case number 
	f = open(gciDir + 'case', 'r')
	case = int(f.readline())
	f.close()

	fig1 = plt.figure()
	fig2 = plt.figure()
	fig3 = plt.figure()
	for i in range(len(dist)):
		(y_ex, u_ex, sigma_ex, 
			u, sigma, 
			k_ex, ksigma_ex, 
			k, sigmak, 
			yc, c_ex, sigmac_ex, 
			c, sigmac, 
			M, 
			Mk,
			Mc
			) = handle(
			case = case, 
			expDir = dist[i], 
			gciDir = gciDir,
			subDir = subDirs[i], 
			alpha = alpha, 
			beta = beta)

		# plot Ux 
		ax1 = fig1.add_subplot(np.ceil(len(dist)/2), 2, i+1) 
		ax1.errorbar(y_ex, u, yerr=sigma)
		ax1.errorbar(y_ex, u_ex, yerr=sigma_ex)
		ax1.set_title('M = ' + str(M))

		# plot k 
		ax2 = fig2.add_subplot(np.ceil(len(dist)/2), 2, i+1)
		ax2.errorbar(y_ex, k, yerr=sigmak)
		# ax2.plot(y_ex, k)
		ax2.errorbar(y_ex, k_ex, yerr=ksigma_ex)
		ax2.set_title('M = ' + str(Mk))

		# plot C 
		ax3 = fig3.add_subplot(np.ceil(len(dist)/2), 2, i+1)
		ax3.plot(yc, c)
		plt.plot(yc, c_ex)
		ax3.set_title('M =' + str(Mc))


	plt.show()