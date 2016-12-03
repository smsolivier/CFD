#!/usr/bin/env python3 

import numpy as np 

import os 
import sys 

import shutil 

import Timer 

import readFoam as rf 

from readFoam import readExperiment # gets experimental values 

''' runs ./run 3 times and writes y profile at each experimental distance ''' 

# --- set up three runs --- 
# number of volumes
Nx1 = np.array([15, 15, 15])
Nx2 = np.array([40, 30, 20])
Ny = np.array([30, 20, 10])
Nz = np.array([30, 20, 10]) 

# total volumes for each run 
N = (Nx1+Nx2)*2*Ny*Nz # number of volumes 

# --- select case --- 
''' case
	0: N339, .6 m/s, same density 
	1: N337, 1 m/s, same density 
	2: N320, .6 m/s, different density 
''' 
# select case, use command line input if provided 
case = 0
if (len(sys.argv) == 2):
	case = sys.argv[1] 

# processor decomposition 
npx = 4 # processors in x 
npy = 1 # processors in y 
npz = 1 # processors in z 

tt = Timer.timer() # start timer 

def handleOF(xinterp, grid):
	''' read OpenFoam output with readFoam at xinterp 
		interpolate onto the specified grid 
	''' 

	y, u, k, C = rf.parse(xinterp)

	uinterp = np.interp(grid, y, u) 

	return uinterp 

def readOF(xinterp):
	''' read openfoam data to extract Ux profile at xinterp ''' 
	y, u, k, C = rf.parse(xinterp)

	return y, u, k, C 

d = np.array([.05, .15, .25, .35, .45]) # distances to interpolate 

# make directory gciData 
dataDir = 'gciData'
if (not(os.path.isdir(dataDir))):
	os.makedirs(dataDir)

# make sub directories for each length 
for i in range(len(d)):
	currentDir = dataDir + '/' + str(d[i])[2:] 
	if (os.path.isdir(currentDir)):
		shutil.rmtree(currentDir) 
	os.makedirs(currentDir) 

# write case to file 
f = open(dataDir + '/case', 'w')
f.write(str(case))
f.close()

# loop backwards so largest run is last (can run paraview on best run) 
for i in range(len(N)-1, -1, -1): # loop through three runs 
	runstr = './run -N {} {} {} {} -quiet -case {} -np {} {} {}'.format(
		int(Nx1[i]), # inlet x 
		int(Nx2[i]), # mixing x 
		int(Ny[i]), # y vols in each half 
		int(Nz[i]), # total z vols 
		case, # which case to run 
		npx, npy, npz # processor decomposition 
		)

	# run openfoam 
	x = os.system(runstr)
	if (x != 0): # exit if problems 
		sys.exit()

	# print all data points at each distance 
	for j in range(len(d)):
		currentDir = dataDir + '/' + str(d[j])[2:] 

		y, u, k, C = readOF(d[j]) 
		# write to file 
		f = open(currentDir+'/run' + str(i) + '.txt', 'w')
		f.write(str(N[i]) + '\n') # write number of volumes at top of file 
		for l in range(len(y)):
			# write y, u, k, c 
			f.write('{:<10e} {:<10e} {:<10e} {:<10e} \n'.format(
				y[l], # grid points 
				u[l], # velocities 
				k[l], # turbelent kinetic energy 
				C[l] # concentration 
				)
			)
		f.close()

tt.stop() # stop timer 