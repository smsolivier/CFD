#!/usr/bin/env python3 

import numpy as np 
import matplotlib.pyplot as plt 

import os 
import sys 

import readFoam as rf 

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
for i in range(len(dist)):
	# read in experimental data 
	df = np.loadtxt(pre+dist[i]+post, skiprows=5, usecols=(0,1,2)) # x, y, Vx

	if (i==0):
		y_ex = np.zeros((len(dist), np.shape(df)[0]))
		u_ex = np.zeros((len(dist), np.shape(df)[0]))

	y_ex[i,:] = df[:,1]
	u_ex[i,:] = -df[:,2] 

def calcGCI(f1, f2, f3, r12, r23, tol=1e-6):
	converged = 0 
	pold = 2 
	alpha = np.fabs((f3 - f2)/(f2 - f1)) 
	while (not(converged)):
		p = np.log(alpha* (r12**pold - 1)/(r23**pold - 1))/np.log(r12)

		if (np.fabs(p - pold)/p < tol):
			converged = 1 

		pold = p 

	return p 

f = open('gci.txt', 'r')
r = f.readline().split() # refinement factors 
f.close()
r12 = float(r[0])
r23 = float(r[1])

ef = np.loadtxt('gci.txt', skiprows=1)

f1 = np.mean(ef[:,1]) 
f2 = np.mean(ef[:,2])
f3 = np.mean(ef[:,3]) 

p = calcGCI(f1, f2, f3, r12, r23)

Fs = 3 
if (np.fabs(p - 2)/2 < .1):
	Fs = 1.25 

sigma = Fs/(r12**p - 1)*np.fabs(f2 - f1)

print(p, sigma)
