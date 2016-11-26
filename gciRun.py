#!/usr/bin/env python3 

import numpy as np 
import matplotlib.pyplot as plt 

import os 
import sys 

import Timer 

import readFoam as rf 

tt = Timer.timer()

def handleOF(xinterp, grid):

	y, u = rf.parse(xinterp)

	uinterp = np.interp(grid, y, u) 

	return uinterp 

grid = np.linspace(-.025, .025, 20)

x1 = np.array([10, 10, 10])
x2 = np.array([40, 30, 20])
y = np.array([40, 30, 20])
z = np.ones(3)*1 

N = (x1+x2)*y*z # number of volumes 

r12 = (N[0]/N[1])**(1/2)
r23 = (N[1]/N[2])**(1/2)

u_of = np.zeros((len(N), len(grid)))
for i in range(len(N)):
	x = os.system('./run -N %s %s %s %s -quiet' % (int(x1[i]), int(x2[i]), int(y[i]), int(z[i])))
	if (x != 0):
		sys.exit()

	u_of[i,:] = handleOF(.05, grid)

# print to file 
f = open('misc/gci.txt', 'w')
f.write(str(r12) + ' ' + str(r23) + '\n')
fmt = '{:<10e}'
for i in range(len(grid)):
	f.write(fmt.format(grid[i]) + ' ') 

	for j in range(3):
		f.write(fmt.format(u_of[j,i]) + ' ')

	f.write('\n')

f.close()

tt.stop()