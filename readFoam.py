#!/usr/bin/env python3 

import numpy as np 
import matplotlib.pyplot as plt 

import os 
import sys 

def parse(xinterp, PLOT=False):
	''' parses OpenFOAM time directories and interpolates the y 
		velocity profile using the writeCellCentres output. 

		Inputs:
			xinterp:
				x value to interpolate at 
			PLOT: 
				plots the Ux profile if true 

		Outputs:
			y1: 
				y values where velocity is calculated 
			uinterp:
				interpolated Ux points 
	''' 

	tdir = '0'
	# find highest time directory 
	for f in os.listdir('.'):
		if (f[0].isdigit() and os.path.isdir(f)):
			if (float(f) > float(tdir)):
				tdir = f 
	tdir += '/'

	print('using time directory ' + tdir)

	initDir = '0/' # place where cell centers are stored 

	ufile = open(tdir+'U', 'r') # open output velocity file 


	# --- get velocities --- 
	# store velocity components into U 
	for line in ufile:
		if (line.startswith('internalField')):
			N = int(next(ufile).strip()) # get number of volumes 
			# skip line 
			next(ufile)

			U = np.zeros((N,3)) # store Ux, Uy, Uz 
			i = 0 # store line number 
			for line in ufile: 
				if (line.strip() == ')'): # break if end of internalField output 
					break 

				# get velocity vector 
				U[i,:] = line[1:-2].split() # remove ( and ), split three components 

				i += 1 # update array location 

	ufile.close() # close velocity file 


	# --- get position --- 
	files = ['ccx', 'ccy', 'ccz'] # name of cell center files in 0 dir 

	r = np.zeros((N,len(files))) # stores position vector 

	for i in range(len(files)): # loop through cell center files 
		g = open(initDir+files[i], 'r') # open file i 
		for line in g: # loop through 
			if (line.startswith('internalField')): # skip until internalField 
				n = int(next(g).strip()) # get number of volumes 
				next(g) # skip line 

				assert (n == N) # make sure sizes same 

				j = 0 # array storage location 
				for line in g:
					if (line.strip() == ')'): # break if end of internalField 
						break 

					# get coordinate 
					r[j,i] = float(line.strip()) # store location 

					j += 1 # update array location 
		g.close()


	# --- interpolate to xinterp --- 
	# get two nearest points, want them to straddle xinterp 
	indx = np.argmin(np.fabs(r[:,0] - xinterp))
	x1 = r[indx,0] # closest value 
	x2 = r[indx+1,0] # next value up 
	if (x1 > xinterp): # make sure x1 and x2 straddle xinterp 
		x2 = r[indx-1,0]

	print('interpolating U at', xinterp, 'between', x1, x2)

	# get z index for 3D 
	zwidth = .025 # width of domain in z 
	indz = np.argmin(np.fabs(r[:2] - zwidth))
	z = r[indz, 2] # closest z value to zwidth 

	y1 = [] # store values at x1 
	y2 = [] # store values at x2

	u1 = [] # store Ux values at x = x1 
	u2 = [] # store Ux values at x = x2 

	for i in range(N):
		# test for x = x1 
		if (r[i,0] == x1 and r[i,2] == z):
			y1.append(r[i,1]) # append y value 
			u1.append(U[i,0]) # append Ux 

		# test for x = x2 
		elif (r[i,0] == x2 and r[i,2] == z):
			y2.append(r[i,1]) # store y value 
			u2.append(U[i,0]) # store Ux 

	# convert to numpy arrays 
	y1 = np.array(y1)
	y2 = np.array(y2)

	u1 = np.array(u1)
	u2 = np.array(u2) 

	# interpolate 
	# y = (y2 - y1)/(x2 - x1) * (x - x1) + y1 
	uinterp = (u2 - u1)/(x2 - x1) * (xinterp - x1) + u1 

	if (PLOT):
		plt.plot(y1*1000, uinterp)
		plt.xlabel('y (mm)')
		plt.ylabel('Ux (m/s)')
		plt.title(str(xinterp))
		plt.show()

	return y1, uinterp 

if __name__ == '__main__':

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
	norm = np.zeros(len(dist))
	plt.figure(figsize=(16,12))
	for i in range(len(dist)):
		# read in experimental data 
		df = np.loadtxt(pre+dist[i]+post, skiprows=5, usecols=(0,1,2)) # x, y, Vx

		y_ex = df[:,1]
		u_ex = -df[:,2] 

		y, u = parse(d[i])

		uinterp = np.interp(y_ex, y*1000, u) 

		norm[i] = np.linalg.norm(uinterp - u_ex, 2)

		plt.subplot(np.ceil(len(dist)/2), 2, i+1)
		plt.plot(y_ex, u_ex, '-o')
		plt.plot(y*1000, u)
		# plt.plot(y_ex, uinterp, '-o')

		plt.title(dist[i])

	for i in range(len(dist)):
		print(dist[i], norm[i])

	plt.show()



