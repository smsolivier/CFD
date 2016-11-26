#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import sys 

import time 

import warnings 
warnings.filterwarnings('ignore') # turn off annoying warnings for mpl 

''' Plots residuals from OpenFOAM stdout ''' 

def getTime():
	try:

		def skipHeader(fname):
			f = open(fname, 'r')
			line = next(f)
			while ('Starting time loop' not in line):
				line = next(f)

			return f 

		fname = 'log'

		var = ['deltaT'] # store variable names 

		# get variable names 
		# f = skipHeader(fname)
		# for line in f:
		# 	if (line.startswith('Time = 1')):
		# 		line = next(f) 
		# 		while ('ExecutionTime' not in line):
		# 			for i in range(len(line)):
		# 				string = 'Solving for '
		# 				if (line[i:i+len(string)] == string):
		# 					var.append(line[i+len(string):].split(',')[0])
		# 			line = next(f)
		# 		break 
		# f.close()

		num = len(var)

		fig = plt.figure(figsize=(8,6))
		# plt.ion()
		run = 1
		while(run == 1):

			val = [] # store values 
			t = [] # store time 
			f = skipHeader(fname)
			for line in f:
				if ('deltaT = ' in line):
					for j in range(len(line)):
						string = 'deltaT = ' 
						if (line[j:j+len(string)] == string):
							val.append(line[j+len(string):].split(',')[0])
				elif (line.startswith('Time = ')):
					l = line.split() 
					t.append(float(l[2]))
				if (line == 'End\n'):
					run = 0

			# print(str(t[-1]) + '\t' + str(val[-1].strip()) + '\r', end='', flush=True)
			# time.sleep(1)

			plt.semilogy(t, val)
			plt.xlabel('t')
			plt.ylabel(r'$\Delta$t')
			plt.title(str(t[-1]))
			f.close()

			plt.legend(loc='best', frameon=False)

			if (run == 1):
				plt.pause(1)
				fig.clear()
			if (run == 0):
				plt.show()

	# exit nicely on ctrl-c 
	except KeyboardInterrupt: 
		sys.exit() 

if __name__ == '__main__':
	getTime()