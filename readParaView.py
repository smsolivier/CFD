#!/usr/bin/env python3 

import numpy as np 
import matplotlib.pyplot as plt 

def read(fname, *args):

	f = open(fname, 'r') # open file 

	header = f.readline().split(',') # read first line 

	f.close() # close file 

	ind = [] # store column numbers of interest 
	names = [] # store names of columns of interest 

	# remove quotes from variable names 
	for i in range(len(header)):
		header[i] = header[i].replace('"', '')

	for i in range(len(args)):
		for j in range(len(header)):
			val = header[j]

			if (val == args[i]):
				ind.append(j)
				names.append(val)

	# load columns of interest 
	df = np.loadtxt(fname, skiprows=1, delimiter=',')

	newdf = np.zeros((np.shape(df)[0], len(ind)))

	for i in range(len(ind)):
		newdf[:,i] = df[:,ind[i]]

	return newdf, names 


if __name__ == '__main__':
	df, names = read('050.csv', 'U:0', 'Points:1')

	print(names)

	plt.plot(df[:,1], df[:,0])
	plt.show()