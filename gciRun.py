#!/usr/bin/env python3 

import numpy as np 
import matplotlib.pyplot as plt 

import os 
import sys 

import shutil 

import Timer 

import readFoam as rf 

tt = Timer.timer()

def handleOF(xinterp, grid):
	''' read OpenFoam output with readFoam at xinterp 
		interpolate onto the specified grid 
	''' 

	y, u = rf.parse(xinterp)

	uinterp = np.interp(grid, y, u) 

	return uinterp 

def readOF(xinterp):
	y, u = rf.parse(xinterp)

	return y, u 

# common interpolation grid to compare to 
grid = np.linspace(-.025, .025, 30)

d = np.array([.05, .15, .25, .35, .45]) # distances to interpolate 

# set up three runs 
Nx1 = np.array([15, 15, 15])
Nx2 = np.array([60, 50, 40])
Ny = np.array([40, 30, 20])
# Nx2 = np.array([10, 10, 10])
# Ny = np.array([10,10,10])
Nz = np.ones(3)*1 

# total volumes for each run 
N = (Nx1+Nx2)*2*Ny*Nz # number of volumes 

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

u_of = np.zeros((len(N), len(d), len(grid))) # store interpolated velocities 
for i in range(len(N)): # loop through three runs 
	x = os.system('./run -N %s %s %s %s -quiet' % (int(Nx1[i]), int(Nx2[i]), int(Ny[i]), int(Nz[i])))
	if (x != 0): # exit if problems 
		sys.exit()

	for j in range(len(d)):
		currentDir = dataDir + '/' + str(d[j])[2:] 

		y, u = readOF(d[j]) 
		# write to file 
		f = open(currentDir+'/run' + str(i) + '.txt', 'w')
		f.write(str(N[i]) + '\n') # write number of volumes at top of file 
		for k in range(len(y)):
			f.write('{:<10e} {:<10e} \n'.format(y[k], u[k]))
		f.close()

		u_of[i,j,:] = handleOF(d[j], grid) 

for i in range(len(d)):

	currentDir = dataDir + '/' + str(d[i])[2:] 

	f = open(currentDir+'/gci.txt', 'w')
	f.write(str(N[0]) + ' ' + str(N[1]) + ' ' + str(N[2]) + '\n') # write number of volumes at top 
	fmt = '{:<10e}'
	for j in range(len(grid)):
		f.write(fmt.format(grid[j]) + ' ') 

		for k in range(3): # loop through 3 runs 
			f.write(fmt.format(u_of[k,i,j]) + ' ')

		f.write('\n')

	f.close()

tt.stop()