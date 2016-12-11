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
# number of volumes, most volumes to least 
Nx1 = np.array([20, 20, 20])
Nx2 = np.array([60, 40, 30])
Ny = np.array([20, 20, 20])
Nz = np.array([1, 1, 1])

# total volumes for each run 
N = (Nx1+Nx2)*2*Ny*Nz # number of volumes 

# --- select case --- 
''' case
	0: N339, .6 m/s, same density 
	1: N337, 1 m/s, same density 
	2: N320, .6 m/s, different density 
''' 
# select case, use command line input if provided 
case = 2
if (len(sys.argv) == 2):
	case = sys.argv[1] 

# processor decomposition 
npx = 4 # processors in x 
npy = 1 # processors in y 
npz = 1 # processors in z 

# set mass diffusivity 
Dab = 1e-6

tt = Timer.timer() # start timer 

# def handleOF(xinterp, grid):
# 	''' read OpenFoam output with readFoam at xinterp 
# 		interpolate onto the specified grid 
# 	''' 

# 	y, u, k, C = rf.parse(xinterp)

# 	uinterp = np.interp(grid, y, u) 

# 	return uinterp 

# def readOF(xinterp):
# 	''' read openfoam data to extract Ux profile at xinterp ''' 
# 	y, u, k, C = rf.parse(xinterp)

# 	return y, u, k, C 

d = np.array([.05, .15, .25, .35, .45]) # distances to interpolate 

# make directory gciData 
dataDir = 'gciData'
if (os.path.isdir(dataDir)):
	shutil.rmtree(dataDir)
os.makedirs(dataDir)

# write case to file 
f = open(dataDir + '/case', 'w')
f.write(str(case))
f.close()

# loop backwards so largest run is last (can run paraview on best run) 
for i in range(len(N)-1, -1, -1): # loop through three runs 
	currentDir = dataDir + '/' + str(N[i]) + '/' # current directory of run 
	# check if exists
	if (os.path.isdir(currentDir)):
		shutil.rmtree(currentDir) # overwrite if it exists 
	os.makedirs(currentDir) # remake 

	runstr = './run \
		-N {} {} {} {} \
		-quiet \
		-case {} \
		-np {} {} {} \
		-Dab {} \
		-outDir {} >gciLog'.format(
		int(Nx1[i]), # inlet x 
		int(Nx2[i]), # mixing x 
		int(Ny[i]), # y vols in each half 
		int(Nz[i]), # total z vols 
		case, # which case to run 
		npx, npy, npz, # processor decomposition 
		Dab, # mass diffusivity 
		currentDir # output directory 
		)

	# run openfoam 
	x = os.system(runstr)
	if (x != 0): # exit if problems 
		sys.exit()

	# print volume distribution to file 
	f = open(currentDir + 'nvols', 'w')
	f.write('{} {} {} {}\n'.format(Nx1[i], Nx2[i], Ny[i], Nz[i]))
	f.close()

	print('completed run', i)

tt.stop() # stop timer 