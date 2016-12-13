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

import texTools as tex # latex helper functions 

from hidespines import *

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

	# integrate omega(x) dx 
	Omega = integrate.quad(omega, -np.inf, np.inf)

	# plt.figure()
	# plt.axvline(u_ex)
	# plt.axvline(u_sim)
	# plt.plot(x, f_ex(x), label='ex')
	# plt.plot(x, f_sim(x), label='sim')
	# plt.plot(x, omega(x), label='omega')
	# plt.title(str(Omega))
	# plt.legend(loc='best')
	# plt.show()

	# print(u_sim/u_ex, sigma_sim/sigma_ex)
	# Omega, err = quadrature(omega, u_sim-1, u_sim+1)

	# print('Omega =', Omega[0], 'Error =', Omega[1])

	return Omega[0] 
	# return 0

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
	r12 = (N[0]/N[1])**(1/2)
	r23 = (N[1]/N[2])**(1/2) 

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

def readGCI(gciDir, subDir, N, grid, cgrid):
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

	U = np.zeros((3,len(grid))) # store velocities on grid 
	K = np.zeros((3, len(grid))) # store k for each run on each grid point 
	C = np.zeros((3, len(cgrid))) # store C for each run on each grid point 

	# store metrics of interest: centerline Ux, centerline k, integral k, integral Ux
	metricNames = ['U', 'k', 'kint', 'Uint']
	metrics = np.zeros((len(metricNames), 3)) 
	for i in range(3): # loop through number of runs 
		curDir = gciDir + str(N[i]) + '/' + subDir + '/' # directory of the run 
		# read in data 
		# U
		yu, u = np.loadtxt(curDir+'U', unpack=True, usecols=(0,1)) # y, Ux 
		# k 
		yk, k = np.loadtxt(curDir+'k', unpack=True)
		# C
		yc, c = np.loadtxt(curDir+'C', unpack=True)

		# interpolate onto grid 
		Ufunc = interp1d(yu, u, kind='cubic') # interpolated function 
		U[i,:] = Ufunc(grid) # evaluated on grid 
		kfunc = interp1d(yk, k, kind='cubic') # interpolated function, used for integration 
		K[i,:] = kfunc(grid) # interpolate k onto grid 
		Cfunc = interp1d(yc, c, kind='cubic') # concentration interpolated function 
		C[i,:] = Cfunc(cgrid) # interpolate onto concentration grid 

		# compute GCI metrics 
		metNum = 0 # which metric to store into 
		metrics[metNum,i] = Ufunc(0) # interpolated centerline Ux 
		metNum += 1
		metrics[metNum,i] = kfunc(0) # interpolated centerline k 
		metNum += 1
		# metrics[2,i] = Cfunc(0) # interpolated centerline concentration 
		# metNum += 1
		metrics[metNum,i] = integrate.quad(kfunc, grid[0], grid[-1])[0] # integral of k 
		metNum += 1
		metrics[metNum,i] = integrate.quad(Ufunc, grid[0], grid[-1])[0] # integral of U 
		metNum += 1

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
		# else: # if a bad metric 
			# print(metricNames[i] + ': ' + bad) # print reason 

	if (len(relError) == 0): # all metrics bad 
		print(subDir, 'no good metrics found')

		error = 0 # set error so that plots ok 

	else: # if at least one good metric 

		relError = np.array(relError) # switch to numpy 
		ind = np.argmax(relError) # location of max error

		error = 2*relError[ind] # maximum relative error, return 2 sigma 

		print(subDir, useMetric[ind], useP[ind], error)

		# make latex table 
		table = tex.table()
		table.addLine(
			subDir, # x location 
			useMetric[ind], # which metric is used 
			tex.utils.writeNumber(useP[ind]), # resulting convergence 
			tex.utils.writeNumber(error) # relative error 
			)
		table.save('report/p_'+subDir+'.tex')

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
	grid = y_ex 

	# create experimental data 
	# interpolate exact answer onto grid 
	ex = u_ex 

	# create simulation data 
	# interpolate simulation onto grid 
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

	# compute error in derivatives 
	E = np.zeros(Ndata)
	for i in range(Ndata): # loop through derivatives 
		if (dsim[i] == 0 and dex[i] == 0):
			ratio = 1 
		elif (dex[i] == 0):
			ratio = 1
		else:
			ratio = dsim[i]/dex[i]

		E[i] = np.fabs(1 - ratio)


	# compute M 
	M = 0
	for i in range(1, len(grid)):
		M += alpha*(1 - Omega[i]) + beta*E[i-1]
	M /= Ndata 

	return M 

def modelError(grid, val_ex, sigma_ex, val, sigma, uinput):
	# interpolate to centerline 
	val_ex_0 = interp1d(grid, val_ex, kind='cubic')(0)
	sigma_ex_0 = interp1d(grid, sigma_ex, kind='cubic')(0)

	val_0 = interp1d(grid, val, kind='cubic')(0) 
	sigma_0 = interp1d(grid, sigma, kind='cubic')(0)

	E = (val_0 - val_ex_0)/val_ex_0

	uval = np.sqrt(sigma_0**2 + sigma_ex_0**2 + uinput**2) 

	return E, uval

def inputUncert(case):
	''' read in input uncertainty ''' 
	caseName = ['N339', 'N337', 'N320','N318'] 

	dr = 'InputUncertainty/' + caseName[case] + '.csv'

	f = open(dr, 'r')

	u = 0 # store overall uncertainty 
	for line in f:
		if line.startswith('Overall'):

			u = float(line.split(',')[1])

	f.close()

	return 2*u # return 2 sigma value 



def handle(case, expDir, gciDir, subDir, N, alpha, beta):
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
	u, sigma, k, sigmak, c, sigmac = readGCI(gciDir, subDir, N, y_ex/1000, yc/1000) 

	# get input uncertainty 
	uinput = inputUncert(case)

	# add in uinput to gci uncertainty 
	combine = lambda x, y: np.sqrt(x**2 + y**2) 
	sigma = combine(sigma, uinput)
	sigmak = combine(sigmak, uinput)
	sigmac = combine(sigmac, uinput)

	# compute M 
	M = globalMerit(y_ex, u_ex, sigma_ex, u, sigma, alpha, beta)

	# compute M for k 
	Mk = globalMerit(y_ex, k_ex, ksigma_ex, k, sigmak, alpha, beta)

	# compute M for C 
	Mc = globalMerit(yc, c_ex, sigmac_ex, c, sigmac, alpha, beta)

	# --- compute model error --- 
	# get relative errors and relative uncertainties 
	Eu, uval_u = modelError(y_ex, u_ex, sigma_ex, u, sigma, uinput)
	Ek, uval_k = modelError(y_ex, k_ex, ksigma_ex, k, sigmak, uinput)
	Ec, uval_c = modelError(yc, c_ex, sigmac_ex, c, sigmac, uinput)

	# L1 norm of three centerline values 
	E = np.linalg.norm([Eu, Ek, Ec], 1) # total relative error 
	uval = np.sqrt(uval_u**2 + uval_k**2 + uval_c**2) # total relative uncertainty 

	return (
		y_ex, u_ex, sigma_ex, # experimental Ux 
		u, sigma, # simulated Ux 
		k_ex, ksigma_ex, # experimental k 
		k, sigmak, # simulated k 
		yc, c_ex, sigmac_ex, # experimental concentration 
		c, sigmac, # simulated concentration 
		M, # metric array for Ux 
		Mk, # metric array for k 
		Mc, # metric for concentration 
		E, # comparison error 
		uval # validation uncertainty 
		)  

if __name__ == '__main__':
	# get command line arguments 
	# undisclosed weighting factors for M formula 
	alpha = 1 
	beta = 1 

	gciDir = 'gciData/'

	if (len(sys.argv) == 2):
		gciDir = sys.argv[1] + '/'

	print('loading', gciDir)
	dist = ['050', '150', '250', '350', '450'] 
	# dist = ['050']
	subDirs = sorted(os.listdir(gciDir)) # get list of subdirectories in gciDir 

	# get case number 
	f = open(gciDir + 'case', 'r')
	case = int(f.readline())
	f.close()

	# get number of volumes on each run 
	ndirs = os.listdir(gciDir) # store directory names inside gciDir 
	N = [] # store volume numbers 
	for i in range(len(ndirs)):
		if (ndirs[i][0].isdigit()):
			N.append(int(ndirs[i]))

	N = np.array(sorted(N, reverse=True)) # sort N from high to low 
	N = N[:3] # only keep largest three

	print(N)

	if (case == 3): # if blind case 
		grid = np.linspace(-.025, .025, 30)
		fig1 = plt.figure()
		fig2 = plt.figure()
		fig3 = plt.figure()
		for i in range(len(dist)):
			u, sigma, k, sigmak, c, sigmac = readGCI(gciDir, dist[i], N, grid, grid)

			# plot Ux 
			ax1 = fig1.add_subplot(len(dist), 1, i+1) 
			ax1.errorbar(grid, u, yerr=sigma, label='simulation')
			plt.xlabel('y (m)')
			plt.ylabel('./g')
			hidespines(ax1)

			# plot k 
			ax2 = fig2.add_subplot(len(dist), 1, i+1)
			ax2.errorbar(grid, k, yerr=sigmak, label='simulation')
			hidespines(ax2)

			# plot C 
			ax3 = fig3.add_subplot(len(dist), 1, i+1)
			ax3.errorbar(grid, c, yerr=sigmac, label='simulation')
			hidespines(ax3)

		# write in requested format 
		# interpolate concentration onto same grid 
		f = open(gciDir+'NorrisOlivier'+dist[i]+'.dat', 'w')
		f.write('Norris, Olivier\n')
		# write header 
		f.write(
			'{:<5} '
			'{:<9} '
			'{:<12} {:<12} {:<12} '
			'{:<12} {:<12} {:<12} '
			'{:<12} {:<12} {:<12} \n'.format(
				'x',
				'y',
				'U', 'U-DU', 'U+DU',
				'K', 'K-DK', 'K+DK',
				'C', 'C-DC', 'C+DC'
				)
			)
		# write values 
		for j in range(len(grid)):
			f.write('{: <5g} ' # x position 
				'{: < 8g} ' # y 
				'{: < 12g} {: < 12g} {: < 12g} ' # velocity 
				'{: < 12g} {: < 12g} {: < 12g} ' # k 
				'{: < 12g} {: < 12g} {: < 12g} \n' # concentration 
				.format(
					float(dist[i])/1000, 
					grid[j], 
					u[j], u[j]-sigma[j], u[j]+sigma[j],
					k[j], k[j] - sigmak[j], k[j] + sigmak[j],
					c[j], c[j]-sigmac[j], c[j]+sigmac[j]
					)
				)
		f.close()

		fig1.savefig('report/blindU.pdf')
		fig2.savefig('report/blindk.pdf')
		fig3.savefig('report/blindC.pdf')
		plt.show()

	else:

		# store M values 
		Mu = np.zeros(len(dist))
		Mk = np.zeros(len(dist))
		Mc = np.zeros(len(dist))

		# store comparison error and uval at each location 
		E = np.zeros(len(dist))
		uval = np.zeros(len(dist))

		# fig1 = plt.figure()
		# fig2 = plt.figure()
		# fig3 = plt.figure()
		for i in range(len(dist)):
			(y_ex, u_ex, sigma_ex, 
				u, sigma, 
				k_ex, ksigma_ex, 
				k, sigmak, 
				yc, c_ex, sigmac_ex, 
				c, sigmac, 
				Mu[i], 
				Mk[i],
				Mc[i],
				E[i],
				uval[i]
				) = handle(
				case = case, 
				expDir = dist[i], 
				gciDir = gciDir,
				subDir = dist[i],
				N = N, 
				alpha = alpha, 
				beta = beta
				)

			# # plot Ux 
			# ax1 = fig1.add_subplot(np.ceil(len(dist)/2), 2, i+1) 
			# ax1.errorbar(y_ex, u, yerr=sigma, label='simulation')
			# ax1.errorbar(y_ex, u_ex, yerr=sigma_ex, label='experiment')
			# ax1.set_title(dist[i])
			# ax1.legend(loc='best', frameon=False)
			# hidespines(ax1)

			# # plot k 
			# ax2 = fig2.add_subplot(np.ceil(len(dist)/2), 2, i+1)
			# ax2.errorbar(y_ex, k, yerr=sigmak, label='simulation')
			# ax2.errorbar(y_ex, k_ex, yerr=ksigma_ex, label='experiment')
			# ax2.set_title(dist[i])
			# ax2.legend(loc='best', frameon=False)
			# hidespines(ax2)

			# # plot C 
			# ax3 = fig3.add_subplot(np.ceil(len(dist)/2), 2, i+1)
			# ax3.errorbar(yc, c, yerr=sigmac, label='simulation')
			# ax3.errorbar(yc, c_ex, yerr=sigmac_ex, label='experiment')
			# # ax3.set_title('M =' + str(Mc[i]))
			# ax3.set_title(dist[i])
			# ax3.legend(loc='best', frameon=False)
			# hidespines(ax3)

			plt.figure()
			plt.subplot(3,1,1)
			plt.errorbar(y_ex, u, yerr=sigma, label='simulation')
			plt.errorbar(y_ex, u_ex, yerr=sigma_ex, label='experiment')
			plt.xlabel('y (mm)')
			plt.ylabel('U')
			hidespines(plt.gca())

			plt.subplot(3,1,2)
			plt.errorbar(y_ex, k, yerr=sigmak, label='simulation')
			plt.errorbar(y_ex, k_ex, yerr=ksigma_ex, label='experiment')
			plt.xlabel('y (mm)')
			plt.ylabel('k')
			hidespines(plt.gca())

			plt.subplot(3,1,3)
			plt.errorbar(yc, c, yerr=sigmac, label='simulation')
			plt.errorbar(yc, c_ex, yerr=sigmac_ex, label='experiment')
			plt.xlabel('y (mm)')
			plt.ylabel('C')
			plt.xlim(yc[-1], yc[0])
			hidespines(plt.gca())

			plt.savefig('report/three_'+dist[i]+'.pdf')

		# fig1.savefig('report/U.pdf')
		# fig2.savefig('report/k.pdf')
		# fig3.savefig('report/C.pdf')

		# L1 norm of all distances 
		Etotal = np.linalg.norm(E, 1) 

		uval_tot = 0
		for i in range(len(dist)):
			uval_tot += uval[i]**2 
		uval_tot = np.sqrt(uval_tot)

		print('\nM Values:')
		print('{:<6} {:<20} {:<20} {:<20}'.format('Dist', 'Mu', 'Mk', 'Mc'))
		table = tex.table()
		for i in range(len(dist)):
			print('{:<6} {:<20} {:<20} {:<20}'.format(dist[i], Mu[i], Mk[i], Mc[i]))
			table.addLine(
				dist[i], # name of line
				tex.utils.writeNumber(Mu[i]), # Mu 
				tex.utils.writeNumber(Mk[i]), # Mk 
				tex.utils.writeNumber(Mc[i]) # Mc 
				)
		table.save('report/M.tex')

		print('\nComparison Error')
		print('{:6} {:20} {:20} {:20}'.format('Dist', 'E', 'uval', 'uval/E'))
		table = tex.table()
		for i in range(len(dist)):
			print('{:<6.10} {:<20.10e} {:<20.10e} {:<20.10}'.format(
				dist[i], 
				E[i], 
				uval[i], 
				uval[i]/np.fabs(E[i])
				)
			)
			table.addLine(
				dist[i], # name of line 
				tex.utils.writeNumber(E[i]), # comparison error 
				tex.utils.writeNumber(uval[i]), # uncertainty 
				tex.utils.writeNumber(uval[i]/np.fabs(E[i])) # ratio 
				)

		# add line for overal value 
		table.addLine(
			'Overall', # name of line 
			tex.utils.writeNumber(Etotal), # comparison error 
			tex.utils.writeNumber(uval_tot), # uncertainty 
			tex.utils.writeNumber(uval_tot/Etotal) # ratio 
			)
		table.save('report/CE'+str(case)+'.tex')


		macro = tex.macro()
		macro.define('Etot'+str(case), Etotal, '{:.3g}')
		macro.define('uval_tot'+str(case), uval_tot, '{:.3g}')
		macro.define('ratio'+str(case)+'.tex', uval_tot/Etotal, '{:.3g}')
		macro.save('report/Emacro_'+str(case)+'.tex')

		print(Etotal, uval_tot, uval_tot/Etotal)

		# plt.show()