#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

test = 1

if (test == 0):
	pre = 'N320/'
	post = '_PIV_dr1_u06.dat'
elif (test == 1):
	pre = 'N337/'
	post = '_PIV_dr0_u10.dat'
elif (test == 2):
	pre = 'N339/'
	post = '_PIV_dr0_u06.dat'

dist = ['050', '150', '250', '350', '450']
# dist = ['050']

for i in range(len(dist)):
	df = np.loadtxt(pre+'X'+dist[i]+post, skiprows=5, usecols=(0,1,2,4)) # x, y, Vx, Vy 

	print('x =', df[0,0])

	plt.plot(df[:,1], -df[:,2], label=dist[i])

plt.title(pre[:-1])
plt.xlabel('y')
plt.ylabel('Vx')
plt.legend(loc='best')
plt.show()
