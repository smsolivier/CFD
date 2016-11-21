#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
from scipy.interpolate import griddata
import os
import sys
import subprocess

# Make an initial condition to check for the 'writeCellCentres' files
# If it is not there the script should be run
# If it is there, the script should be verified to comply with the
#   'blockMeshDict' file, this could check for the total number of cells
#   at the face.

# Read in the data in a useable format
def readin(filename, skip=0):
	"""Reads in the output file"""
	f = open(filename, 'r')
	with open(filename) as f:
	    lines = f.readlines()
	f.close()
	return lines[skip:]

def grabValues(file, str):
	vals = []
	# Go to the line immediately after 'str'
	for i, line in enumerate(file):
		if (line.strip() == str):
			break
	# This has is now at the line immediately after 'str'
	# The 'writeCellCentres' can be written in two forms:
	#   The first is for few cells (~10), everything is written on the same line
	#   The second is for many cells, everything is writtien in a list
	# The form can be determined by checking the final char of the 'values...' lines
	# If the char is ';' it is of the FIRST FORM, if it is '>' it is of the SECOND FORM
	FORM = 1 # Assumes the file is of the first form
	for line in file:
		# The first if-statement searches for the first instance of 'value' following the
		# patches keyword
		if (line.strip()[:5] == 'value'):
			# The next if-statements differentitate between the two form
			# This is for the FISRT FORM
			if (line.strip()[-1] == ';'):
				# This line must be immediately parsed to determine the values
				# before the outer for-loop passes over the line
				for i, char in enumerate(line.strip()):
					# This for-loop looks through each character of the line and
					# determines the indicies between where the values are stored
					if (char == '('):
						beginning = i+1
					elif (char == ')'):
						end = i
					else:
						continue
				# With the index values determine, the values are extracted into a list
				vals = np.fromstring(line.strip()[beginning:end], sep=' ')
				# With the values found, the loop can be exited and values returned
				break
			# This is for the SECOND FORM
			elif (line.strip()[-1] == '>'):
				# The current line is not necessary for this form
				# So the loop is exited and continued below
				FORM = 2
				break
			# This throws an error if the form cannot be determined...
			else:
				print("Error. Could not determine which form the cellCentres are saved as... ask Jesse.")
				break

	if (FORM == 2):
		read = False
		for line in file:
			# Continuation of the SECOND FORM loop
			# Reads through every line following the 'value...' line
			# Saves all the values between the parenthesis to the 'vals' list
			# It is of Form 2
			if (line == ')\n'):
				# Exits the loop once it finds the close parenthesis
				break
			if (read == True):
				vals.append(float(line))
			if (line == '(\n'):
				# Begins the loop once it finds the open parenthesis
				read = True
	file.seek(0) # Resets the read position of the file
	return vals

factorV = 10
factorK = 1
PLOT = True
LOUD = False
conversion = 1000
### Project file directory
projectPATH = os.path.dirname(os.path.realpath(__file__))
# projectPATH = '/home/jesse/OpenFOAM/jesse-4.1/run/Project/'
cellCenters = np.array(['ccx', 'ccy', 'ccz'])
### Patch identifiers
patches = np.array(['topInlet', 'bottomInlet'])
### Data path
dataPATH = '/data'
data = '/Vel_Inlet_X-50_u06ms.dat'
### Initial Velocity Path
initialPATH = '/0/'

# Find the location of the cell centered files (either in '/0' or '/constant')
# Removes them once they are found and 
FOUND = False
print("Checking for cell center file locations... ", end='')
for root, dirs, files in os.walk(projectPATH):
	for file in files:
		# print(os.path.join(root, file))
		if (file == cellCenters[0]):
			if ((('constant' in root) or ('0' in root))
					and ('processor' not in root)
					and ('0.' not in root)):
				FOUND = True
				centersPATH = root.replace(projectPATH, '') + '/'
if (FOUND):
	print("Found.")
	print("Removing previous cell center postions... ", end='')
	for coord in cellCenters:
		os.system('rm {}'.format(projectPATH+centersPATH+coord))
	print("Done!")
	print("Rerunning 'writeCellCentres'... ")
else:
	print("Not found.")
	print("\nRunning 'writeCellCentres'...", end='')
os.system('writeCellCentres')
REMOVE = True
print("Done!")
print("Rechecking cell center file locations... ", end='')
for root, dirs, files in os.walk(projectPATH):
	for file in files:
		# print(os.path.join(root, file))
		if (file == cellCenters[0]):
			if ((('constant' in root) or ('0' in root))
					and ('processor' not in root)
					and ('0.' not in root)):
				centersPATH = root.replace(projectPATH, '') + '/'
print("Done!")

listvals = []
# Read the file and gather the coordinates into a list
for coord, name in enumerate(cellCenters):
	print("Reading cell centers from '{}'... ".format(projectPATH+centersPATH+name), end='')
	current = open(projectPATH+centersPATH+name, 'r')
	for pos, patch in enumerate(patches):
		dummy = grabValues(current, patch)
		listvals.append(dummy)
		if (LOUD == True):
			print("\nFrom '{}' patch in the '{}' file:".format(patch, projectPATH+centersPATH+name))
			print("\t", dummy)
	current.close()
	print("Done!")
#####################################################################
# 'listvals' has all the x, y, and z coordinates for the top & bottom inlets
# stored individually, in list form
# Indicies:
#   0 - top inlet x values
#   1 - bottom inlet x values
#   2 - top inlet y values
#   3 - bottom inlet y values
#   4 - top inlet z values
#   5 - bottom inlet z values
#####################################################################

nTopCells = len(listvals[0])
nBottomCells = len(listvals[1])
nCells = nTopCells + nBottomCells

print("Saving cell center positions in a matrix... ", end='')
matvals = np.zeros([nCells, 3])
# Convert the list coordinates into matrix form
for i in range(3):
	# print(np.concatenate(listvals[2*i:2*(i+1)]))
	matvals[:,i] = np.concatenate(listvals[2*i:2*(i+1)])*conversion
print("Done!")
if (LOUD == True):
	print("\tThere are a total of {} cell centers ({} for '{}' and {} for '{}')".format(
		nCells, nTopCells, patches[0], nBottomCells, patches[1]))

#####################################################################
# 'matvals' has all the x, y, and z coordinates for the top & bottom inlets
# stored together, in matrix form
# Indicies:
#   0 - x values
#   1 - y values
#   2 - z values
#####################################################################

if (PLOT == True):
	# Plot the data
	fig = plt.figure()
	ay = fig.add_subplot(111, projection='3d')

	#####################################################################
	### Plots the position data
	# Blue has no velocity
	# Green has x and y components
	# Red has x, y, and z components
	# *For plotting purposes the y and z values are swapped, this doesn't
	#  change any of the data
	#####################################################################
	for i in range(nCells):
		x = matvals[i,0]
		y = matvals[i,1]
		z = matvals[i,2]
		ay.scatter(x, z, y, color='#500000', s=5, lw=0.1)
	#####################################################################

	#####################################################################
	# Overlay the inlets defined in my OpenFOAM input
	#####################################################################
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
	#####################################################################

	ay.view_init(elev=20, azim=-25)
	ay.set_xlabel('X (mm)')
	ay.set_ylabel('Z (mm)')
	ay.set_zlabel('Y (mm)')
	ay.auto_scale_xyz([-55, -45], [-40, 40], [-40, 40])

#####################################################################
# Load in the data for interpolation
#####################################################################
print("Reading inlet data from '{}'... ".format(projectPATH+dataPATH+data), end='')
lines = readin(projectPATH+dataPATH+data, skip=5)
cols = np.array(["x", "y", "z", "Vx", "U_x95", "Vy", "U_y95", "Vz", "U_z95", "RMSVx", "RMSVxUpper", "RMSVxLower", "RMSVy", "RMSVyUpper", "RMSVyLower", "RMSVz", "RMSVzUpper", "RMSVzLower", "k", "kLOWER", "kUPPER"])
nRows = len(lines)
nCols = len(cols)

# Create a matrix of all the values from the file
matdata = np.zeros([nRows, nCols])
for i in range(len(lines)):
	matdata[i,:] = np.fromstring(lines[i], sep='\t')

posdata = np.zeros([nRows, 2])
veldata = np.zeros([nRows, 3])
kdata = np.zeros(nRows)
# Extract the needed position and velocity data
for i in range(nRows):
	# print(np.array([matdata[i,1], matdata[i,2]]))
	posdata[i,:] = np.array([matdata[i,1], matdata[i,2]])
	kdata[i] = np.array([matdata[i,18]])
	for k in range(3):
		veldata[i,k] = matdata[i,2*k+3]
	if (LOUD == True): print(posdata[i], veldata[i])
print("Done!")
#####################################################################

#####################################################################
# pos data contains the non-uniform grid data for the position
# veldata contains the corresponding velocity component values
#####################################################################

# print("position data:\n", posdata, end='\n\n')
# print("x-velocity data:\n", veldata[:,0], end='\n\n\n')
# print("y-cell position:\n", matvals[:,1], end='\n\n')
# print("z-cell position:\n", matvals[:,2], end='\n\n\n')
# print((matvals[:,1], matvals[:,2]))
# print(matvals[:,1:2])

# plt.scatter(matvals[:,2], matvals[:,1])
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(121)
# bx = fig.add_subplot(122)
# ax.scatter(posdata[:,:,2], posdata[:,:,1])
# bx.plot(matvals[:,2], matvals[:,1])
# ax.set_title("Given Data")
# bx.set_title("Data Read from OpenFOAM")
# plt.show()

scheme = 'linear'
# Reformat the data points taken from OpenFOAM to work nicely with interpolate
print("Interpolating U inlet conditions using '{}' scheme... ".format(scheme), end='')
n = len(matvals[:,2])
matpos = np.zeros([n,2])
for i in range(n):
	matpos[i] = np.array([matvals[i,1], matvals[i,2]])

newVx = griddata(posdata,
	veldata[:,0],
	matpos, method=scheme)
newVy = griddata(posdata,
	veldata[:,1],
	matpos, method=scheme)
newVz = griddata(posdata,
	veldata[:,2],
	matpos, method=scheme)
print("Done!")

##########################################################################
# V = np.sqrt(newVx**2 + newVy**2 + newVz**2)
# print(matpos[:,1])
# print(matpos[:,0])
# dz = 0.625
# dy = 1.25
# print("dz =", dz)
# print("dy =", dy)
# print("Average Velocity =", np.sum(V*dz*dy)/(25*50))
# print("Average Velocity =", np.mean(V))
##########################################################################

print("Gathering current U inlet conditions from '{}'... ".format(projectPATH+initialPATH+'U'), end='')
lines = readin(projectPATH+initialPATH+'U')
topFound = False
botFound = False
for i, line in enumerate(lines):
	if (line.strip() == patches[0]):
		topStartPos = i+1 # To account for the parenthesis
		topFound = True
		# print("Found the beginning of '{}' at line {}.".format(patches[0], i+1))
	elif (topFound == True and line.strip() == '}'):
		topEndPos = i
		topFound = False
		# print("Found the ending of '{}' at line {}.".format(patches[0], i))
	elif (line.strip() == patches[1]):
		botStartPos = i+1
		botFound = True
		# print("Found the beginning of '{}' at line {}.".format(patches[1], i+1))
	elif (botFound == True and line.strip() == '}'):
		botEndPos = i
		botFound = False
		# print("Found the ending of '{}' at line {}.".format(patches[1], i))
print("Done!")

# run = np.sort(np.concatenate([np.arange(topStartPos+1),
# 	np.arange(topEndPos,botStartPos+1),
# 	np.arange(botEndPos,len(lines))]))

print("Writing interpolated U inlet conditions to '{}'... ".format(projectPATH+initialPATH+'U'), end='')
# Print everything to a new file
# Print the beginning of the file up till the inlet patch conditions
fmt = '%.15f'
f = open(projectPATH+initialPATH+'U', 'w')
for i in np.arange(topStartPos+1):
	f.write(lines[i])
# Write in the new inlet data
f.write('        type            fixedValue;\n')
f.write('        value           nonuniform List<vector>\n')
f.write('{}\n'.format(nTopCells))
f.write('(\n')
for i in range(nTopCells):
	f.write('({} {} {})\n'.format(fmt%(factorV*newVx[i]), fmt%(factorV*newVy[i]), fmt%(factorV*newVz[i])))
f.write(')\n')
f.write(';\n')
# Continue onto the bottom inlet
for i in np.arange(topEndPos,botStartPos+1):
	f.write(lines[i])
# Write in the new inlet data
f.write('        type            fixedValue;\n')
f.write('        value           nonuniform List<vector>\n')
f.write('{}\n'.format(nTopCells))
f.write('(\n')
for i in range(nTopCells,nCells):
	f.write('({} {} {})\n'.format(fmt%(factorV*newVx[i]), fmt%(factorV*newVy[i]), fmt%(factorV*newVz[i])))
f.write(')\n')
f.write(';\n')
# Finish the rest of the file
for i in np.arange(botEndPos,len(lines)):
	f.write(lines[i])
f.close()
print("Done!")

# print("Backing up previous U inlet conditions to '{}'... ".format(projectPATH+initialPATH+'.Ubackup'), end='')
# f = open(projectPATH+initialPATH+'.Ubackup', 'w')
# for line in lines:
# 	f.write(line)
# f.close()
# print("Done!")


# new = open(U.txt, 'w')
# print(nTopCells)
# print(nBottomCells)
# print(lines)
# print(newVx.shape)









# if (PLOT == True):
# 	fig = plt.figure()
# 	az = plt.add_subplot(111)
# 	az.scatter()




# Reformat the data points taken from OpenFOAM to work nicely with interpolate
print("Interpolating k inlet conditions using '{}' scheme... ".format(scheme), end='')
newK = griddata(posdata,
	kdata[:],
	matpos, method=scheme)
print("Done!")

print("Gathering current k inlet conditions from '{}'... ".format(projectPATH+initialPATH+'k'), end='')
lines = readin(projectPATH+initialPATH+'k')
topFound = False
botFound = False
for i, line in enumerate(lines):
	if (line.strip() == patches[0]):
		topStartPos = i+1 # To account for the parenthesis
		topFound = True
		# print("Found the beginning of '{}' at line {}.".format(patches[0], i+1))
	elif (topFound == True and line.strip() == '}'):
		topEndPos = i
		topFound = False
		# print("Found the ending of '{}' at line {}.".format(patches[0], i))
	elif (line.strip() == patches[1]):
		botStartPos = i+1
		botFound = True
		# print("Found the beginning of '{}' at line {}.".format(patches[1], i+1))
	elif (botFound == True and line.strip() == '}'):
		botEndPos = i
		botFound = False
		# print("Found the ending of '{}' at line {}.".format(patches[1], i))
print("Done!")

print("Writing interpolated k inlet conditions to '{}'... ".format(projectPATH+initialPATH+'k'), end='')
# Print everything to a new file
# Print the beginning of the file up till the inlet patch conditions
fmt = '%.15f'
f = open(projectPATH+initialPATH+'k', 'w')
for i in np.arange(topStartPos+1):
	f.write(lines[i])
# Write in the new inlet data
f.write('        type            fixedValue;\n')
f.write('        value           nonuniform List<scalar>\n')
f.write('{}\n'.format(nTopCells))
f.write('(\n')
for i in range(nTopCells):
	f.write('{}\n'.format(fmt%(factorK*newK[i])))
f.write(')\n')
f.write(';\n')
# Continue onto the bottom inlet
for i in np.arange(topEndPos,botStartPos+1):
	f.write(lines[i])
# Write in the new inlet data
f.write('        type            fixedValue;\n')
f.write('        value           nonuniform List<scalar>\n')
f.write('{}\n'.format(nTopCells))
f.write('(\n')
for i in range(nTopCells,nCells):
	f.write('{}\n'.format(fmt%(factorK*newK[i])))
f.write(')\n')
f.write(';\n')
# Finish the rest of the file
for i in np.arange(botEndPos,len(lines)):
	f.write(lines[i])
f.close()
print("Done!")

if (PLOT == True):
	fig = plt.figure()
	ax = fig.add_subplot(421)
	bx = fig.add_subplot(422)
	cx = fig.add_subplot(423)
	dx = fig.add_subplot(424)
	ex = fig.add_subplot(425)
	fx = fig.add_subplot(426)
	gx = fig.add_subplot(427)
	hx = fig.add_subplot(428)

	ax.scatter(posdata[:,1], posdata[:,0], c=veldata[:,0], cmap='plasma', lw=0)
	bx.scatter(matpos[:,1], matpos[:,0], c=newVx, cmap='plasma', lw=0)
	cx.scatter(posdata[:,1], posdata[:,0], c=veldata[:,1], cmap='plasma', lw=0)
	dx.scatter(matpos[:,1], matpos[:,0], c=newVy, cmap='plasma', lw=0)
	ex.scatter(posdata[:,1], posdata[:,0], c=veldata[:,2], cmap='plasma', lw=0)
	fx.scatter(matpos[:,1], matpos[:,0], c=newVz, cmap='plasma', lw=0)
	gx.scatter(posdata[:,1], posdata[:,0], c=kdata[:], cmap='plasma', lw=0)
	hx.scatter(matpos[:,1], matpos[:,0], c=newK, cmap='plasma', lw=0)

	ax.set_title("Before Interpolation")
	bx.set_title("After Interpolation")
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


if (PLOT == True):
	plt.show()