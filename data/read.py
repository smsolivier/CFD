#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

case = 0

''' read in experimental data for location dist ''' 
field = ['_PIV_', '_LIF_'] # velocity/k v concentration 
if (case == 2):
	pre = 'N320/'
	post = 'dr1_u06.dat'
elif (case == 1):
	pre = 'N337/'
	post = 'dr0_u10.dat'
elif (case == 0):
	pre = 'N339/'
	post = 'dr0_u06.dat'

pre = pre + 'X'

dist = ['050', '150', '250', '350', '450']

for i in range(len(dist)):
	df = np.loadtxt(pre+dist[i]+field[0]+post, skiprows=5, usecols=(0,1,2)) # x, y, Vx

	cf = np.loadtxt(pre+dist[i]+field[1]+post, skiprows=3, usecols=(0,1,2))

	plt.plot(df[:,1], -df[:,2], label=dist[i])

plt.title(pre[:-1])
plt.xlabel('y')
plt.ylabel('Vx')
plt.legend(loc='best')
plt.show()
