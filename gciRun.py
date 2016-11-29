#!/usr/bin/env python3 

import numpy as np 
import matplotlib.pyplot as plt 

import os 
import sys 

import shutil 

import Timer 

import readFoam as rf 

''' runs ./run 3 times and writes y profile at each experimental distance ''' 

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

# common interpolation grid to compare to 
grid = np.linspace(-.025, .025, 30)

d = np.array([.05, .15, .25, .35, .45]) # distances to interpolate 

# set up three runs 
# number of volumes
Nx1 = np.array([15, 15, 15])
Nx2 = np.array([40, 30, 20])
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

# u_of = np.zeros((len(N), len(d), len(grid))) # store interpolated velocities 
for i in range(len(N)): # loop through three runs 
	x = os.system('./run -N %s %s %s %s -quiet' % (int(Nx1[i]), int(Nx2[i]), int(Ny[i]), int(Nz[i])))
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
				y[l], 
				u[l],
				k[l], 
				C[l]
				)
			)
		f.close()

		# u_of[i,j,:] = handleOF(d[j], grid) 

# # print interpolated common grid points for each time 
# # puts all three runs on the same grid in one file 
# for i in range(len(d)):

# 	currentDir = dataDir + '/' + str(d[i])[2:] 

# 	f = open(currentDir+'/gci.txt', 'w')
# 	f.write(str(N[0]) + ' ' + str(N[1]) + ' ' + str(N[2]) + '\n') # write number of volumes at top 
# 	fmt = '{:<10e}'
# 	for j in range(len(grid)):
# 		f.write(fmt.format(grid[j]) + ' ') # write the y points 

# 		for k in range(3): # loop through 3 runs 
# 			f.write(fmt.format(u_of[k,i,j]) + ' ') # write Ux for each run 

# 		f.write('\n')

# 	f.close()

tt.stop() # stop timer 