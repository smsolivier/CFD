#!/usr/bin/env python3 

import numpy as np 
import matplotlib.pyplot as plt 

import os 
import sys 

from scipy.interpolate import interp1d 
from scipy.integrate import quadrature 
import scipy.integrate as integrate 
from scipy.stats import norm

import readFoam as rf 

''' processes gciRun data 
	calculates uncertainties using GCI 
	calculates Omega 
''' 

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
	f = lambda x: 1/(np.sqrt(2*sigma**2*np.pi))*np.exp(
		-1*(x-mu)**2/(2*sigma**2))

	return f 

def calcOmega(u_sim, sigma_sim, u_ex, sigma_ex):
	a = 1/(2*sigma_sim*np.sqrt(np.pi))
	b = lambda x: ((x - u_ex)/sigma_ex)**2 
	c = lambda x: ((x - u_sim)/sigma_sim)**2 

	omega = lambda x: a * np.exp(-1/4*(b(x) + c(x)))

	# x = np.linspace(u_sim-2, u_sim+2, 100)
	# plt.plot(x, omega(x))
	# plt.title(str(u_sim/u_ex))
	# plt.show()

	x = np.linspace(u_sim - 2, u_sim + 2, 100)
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

# print(calcOmega(.3, .25, .6, .3))


def process():
	dirName = 'adaData/' # name of file to read GCI run data from 

	dirs = sorted(os.listdir(dirName)) # get list of directory names 

	# get correct file names according to case 
	case = 2

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

	dist = ['050', '150', '250', '350', '450'] # experimental distances
	d = np.array([.05, .15, .25, .35, .45])

	# --- read in experimental data --- 
	# y_ex stores all y values for all distances 
	# u_ex stores Ux for all distances 
	# sigma_ex stores U_x95 
	for i in range(len(dist)):
		# read in experimental data 
		df = np.loadtxt(pre+dist[i]+post, skiprows=5, usecols=(1,2,3)) # y, Vx, U_x95

		if (i==0):
			y_ex = np.zeros((len(dist), np.shape(df)[0]))
			u_ex = np.zeros((len(dist), np.shape(df)[0]))
			sigma_ex = np.zeros((len(dist), np.shape(df)[0]))

		y_ex[i,:] = df[:,0]
		u_ex[i,:] = -df[:,1] 
		sigma_ex[i,:] = df[:,2]

	# read in gciRun data, interpolate to common grid, do gci on each grid point 
	for j in range(len(dirs)): # loop through time directories 
		N = np.zeros(3) 

		grid = y_ex[j,:]/1000
		# grid = np.linspace(-.024, .024, 50)
		U = np.zeros((3,len(grid)))
		for i in range(3): # loop through number of runs 
			fname = dirName + dirs[j] + '/run'+str(i)+'.txt'
			f = open(fname, 'r')
			N[i] = float(f.readline().strip())

			y, u = np.loadtxt(fname, skiprows=1, unpack=True)

			# U[i,:] = np.interp(grid, y, u) # interpolate onto grid 
			U[i,:] = interp1d(y, u, kind='cubic')(grid) 

		r12 = (N[0]/N[1])**(1/3)
		r23 = (N[1]/N[2])**(1/3) 

		n = len(grid) # number of grid points 
		p = np.zeros(n) # store convergence for each velocity point 
		sigma = np.zeros(n) # store uncertainty 
		Omega = np.zeros(n)
		for i in range(n):
			p[i], sigma[i] = calcGCI(U[0,i], U[1,i], U[2,i], r12, r23) 

			Omega[i] = calcOmega(U[0,i], 2*np.fabs(sigma[i]), u_ex[j,i], sigma_ex[j,i])
			print(Omega[i])
			# Omega[i] = calcOmega(U[0,i], 10*sigma[i], u_ex[j,i], 10*sigma_ex[j,i])


		# print(Omega)

		# average velocity GCI 
		f1 = np.mean(U[0,:])
		f2 = np.mean(U[1,:])
		f3 = np.mean(U[2,:])

		pmean, sigmaMean = calcGCI(f1, f2, f3, r12, r23)
		u_mean_ex = np.mean(u_ex[j,:])
		u_sigma_ex = np.mean(sigma_ex[j,:])

		# Omega = calcOmega(f1, np.fabs(sigmaMean), u_mean_ex, u_sigma_ex)
			# print('Omega =', Omega)
		# print(d[j], 'p =', pmean, 'Uavg =', f1, '+-', 2*sigmaMean)

		plt.subplot(np.ceil(len(dirs)/2), 2, j+1)
		plt.errorbar(grid*1000, U[0,:], yerr=sigma)
		plt.plot(y_ex[j,:], u_ex[j,:])
		# plt.plot(y_ex[j,:], Omega)
		plt.title(dirs[j])

	plt.show()

process()