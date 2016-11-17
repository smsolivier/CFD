#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

dist = ['050', '150', '250', '350', '450']
pre = 'N339/'
# post = '_PIV_dr1_u06.dat'
# post = '_PIV_dr0_u10.dat'
post = '_PIV_dr0_u06.dat'

for i in range(len(dist)):
	df = np.loadtxt(pre+'X'+dist[i]+post, skiprows=5, usecols=(0,1,2,4))

	plt.plot(df[:,1], df[:,2], label=dist[i])
plt.legend(loc='best')
plt.show()