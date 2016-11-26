#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3

CMAP = 'jet'
plotData = np.array(['matvals', 'posdata', 'veldata', 'matpos',
	'kdata', 'newVx', 'newVy', 'newVz', 'newK'])
for f in plotData:
	globals()[f] = np.loadtxt('inletPlotData/'+f, delimiter=',')

#####################################################################
# Cell Center Position Plot
#####################################################################
fig2 = plt.figure()
ay = fig2.add_subplot(111, projection='3d')

##########
### Plots the position data
# Blue has no velocity
# Green has x and y components
# Red has x, y, and z components
# *For plotting purposes the y and z values are swapped, this doesn't
#  change any of the data
##########
for i in range(len(matvals[:,0])):
	x = matvals[i,0]
	y = matvals[i,1]
	z = matvals[i,2]
	ay.scatter(x, z, y, color='#500000', s=5, lw=0.1)
##########

##########
# Overlay the inlets defined in my OpenFOAM input
##########
toppos = np.array([(-49.98287, 25, 1.3088),
				(-49.32845, 25, 26.30028),
				(-49.32845, -25, 26.30028),
				(-49.98287, -25, 1.30885)])
botpos = np.array([(-49.98287, 25, -1.30885,),
			   (-49.32845, 25, -26.30028,),
			   (-49.32845, -25, -26.30028,),
			   (-49.98287, -25, -1.30885,)])
top = a3.art3d.Poly3DCollection([toppos])
bot = a3.art3d.Poly3DCollection([botpos])
top.set_color('white')
bot.set_color('white')
top.set_alpha(0.5)
bot.set_alpha(0.5)
top.set_edgecolor('k')
bot.set_edgecolor('k')
ay.add_collection3d(top)
ay.add_collection3d(bot)
##########

ay.view_init(elev=20, azim=-25)
ay.set_xlabel('X (mm)')
ay.set_ylabel('Z (mm)')
ay.set_zlabel('Y (mm)')
ay.set_title('Cell Center Positions Overlayed on the Inlets\n(Visual Verification)')
ay.auto_scale_xyz([-55, -45], [-40, 40], [-40, 40])
#####################################################################

#####################################################################
# Interpolated Value Plot
#####################################################################
fig1 = plt.figure()
ax = fig1.add_subplot(421)
bx = fig1.add_subplot(422)
cx = fig1.add_subplot(423)
dx = fig1.add_subplot(424)
ex = fig1.add_subplot(425)
fx = fig1.add_subplot(426)
gx = fig1.add_subplot(427)
hx = fig1.add_subplot(428)

ax.scatter(posdata[:,1], posdata[:,0], c=veldata[:,0], cmap=CMAP, lw=0)
bx.scatter(matpos[:,1], matpos[:,0], c=newVx, cmap=CMAP, lw=0)
cx.scatter(posdata[:,1], posdata[:,0], c=veldata[:,1], cmap=CMAP, lw=0)
dx.scatter(matpos[:,1], matpos[:,0], c=newVy, cmap=CMAP, lw=0)
ex.scatter(posdata[:,1], posdata[:,0], c=veldata[:,2], cmap=CMAP, lw=0)
fx.scatter(matpos[:,1], matpos[:,0], c=newVz, cmap=CMAP, lw=0)
gx.scatter(posdata[:,1], posdata[:,0], c=kdata[:], cmap=CMAP, lw=0)
hx.scatter(matpos[:,1], matpos[:,0], c=newK, cmap=CMAP, lw=0)

ax.set_title("Given Data")
bx.set_title("Interpolated Mesh Data")
ax.set_ylabel(r"$V_x$")
cx.set_ylabel(r"$V_y$")
ex.set_ylabel(r"$V_z$")
gx.set_ylabel(r"$k$")
ax.set_xlim(-30,30)
bx.set_xlim(-30,30)
cx.set_xlim(-30,30)
dx.set_xlim(-30,30)
ex.set_xlim(-30,30)
fx.set_xlim(-30,30)
gx.set_xlim(-30,30)
hx.set_xlim(-30,30)
#####################################################################

# Remove the files
for f in plotData:
	os.system('rm inletPlotData/{}'.format(f))

plt.show()