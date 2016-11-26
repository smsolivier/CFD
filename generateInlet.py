#!/usr/bin/env python3

import os
import sys
import numpy as np
from scipy.interpolate import griddata

# Make an initial condition to check for the 'writeCellCentres' files
# If it is not there the script should be run
# If it is there, the script should be verified to comply with the
#   'blockMeshDict' file, this could check for the total number of cells
#   at the face.

def error(str):
	print('############################################################')
	print('Error Details:')
	print('\t', str)
	print('############################################################')

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
			# Throws an error if the form cannot be determined...
			else:
				error("Could not determine which form the 'cellCentre' files are saved as... ask Jesse.")
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



def generateInlet(data, scheme, PLOT, STATUS):
	data = ''.join(data)
	scheme = ''.join(scheme)
	
	conversion = 1000

	# Reinitializes the data file name
	if (data == '0.6' or data == '.6'):
		data = 'Vel_Inlet_X-50_u06ms.dat'
	elif (data == '1.0' or data == '1'):
		data = 'Vel_Inlet_X-50_u10ms.dat'
	else:
		error('Invalid inlet conditions, either \'0.6\' or \'1.0\'.')

	# Determine the directories
	projectPATH = os.path.dirname(os.path.realpath(__file__))

	cellCenters = np.array(['ccx', 'ccy', 'ccz'])
	### Patch identifiers
	patches = np.array(['topInlet', 'bottomInlet'])

	#####################################################################
	# Find the directory where the 'cellCentre' files are located
	# This is needed because the files are either generated in the '/0/'
	# or the '/constant/' directories... not sure why it changes
	# could depend on the version of OpenFOAM being used.
	#####################################################################
	# Find the location of the cell centered files (either in '/0' or '/constant')
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
					centersPATH = root + '/'

	if (FOUND):
		print("Found.")
		print("Removing previous cell center postions... ", end='')
		for coord in cellCenters:
			os.system('rm {}'.format(centersPATH+coord))
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
			# print(os.path.join(root, file)) # prints every file found
			if (file == cellCenters[0]):
				if ((('constant' in root) or ('0' in root))
						and ('processor' not in root)
						and ('0.' not in root)):
					centersPATH = root + '/'
	print("Done!")
	#####################################################################

	listvals = []
	# Read the file and gather the coordinates into a list
	for coord, name in enumerate(cellCenters):
		print("Reading cell centers from '{}'... ".format(centersPATH+name), end='')
		current = open(centersPATH+name, 'r')
		for pos, patch in enumerate(patches):
			dummy = grabValues(current, patch)
			listvals.append(dummy)
			if (STATUS == True):
				print("\nFrom '{}' patch in the '{}' file:".format(patch, centersPATH+name))
				print("\t", dummy)
				print("\tTotal of {} values.".format(len(dummy)))
		current.close()
		print("Done!\n")

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
		matvals[:,i] = np.concatenate(listvals[2*i:2*(i+1)])*conversion
	print("Done!")
	if (STATUS == True):
		print("There are a total of {} cell centers ({} for '{}' and {} for '{}')".format(
			nCells, nTopCells, patches[0], nBottomCells, patches[1]))

	#####################################################################
	# 'matvals' has all the x, y, and z coordinates for the top & bottom inlets
	# stored together, in matrix form
	# Indicies:
	#   0 - x values
	#   1 - y values
	#   2 - z values
	#####################################################################

	#####################################################################
	# Load in the data for interpolation
	#####################################################################
	print("Reading inlet data from '{}'... ".format(projectPATH+'/data/'+data), end='')
	lines = readin(projectPATH+'/data/'+data, skip=5)
	cols = np.array(["x", "y", "z", "Vx", "U_x95", "Vy", "U_y95", "Vz", "U_z95",
		"RMSVx", "RMSVxUpper", "RMSVxLower", "RMSVy", "RMSVyUpper", "RMSVyLower",
		"RMSVz", "RMSVzUpper", "RMSVzLower", "k", "kLOWER", "kUPPER"])
	# Remove the black lines at the end of the file
	for i in range(len(lines))[::-1]:
		if (lines[i].strip() == ''):
			del lines[-1]
		else:
			break
	nRows = len(lines)
	nCols = len(cols)

	# Create a matrix of all the values from the file
	matdata = np.zeros([nRows, nCols])
	for i in range(len(lines)):
		matdata[i,:] = np.fromstring(lines[i], sep='\t')

	##########
	# 'posdata' contains the non-uniform grid data for the position
	# 'kdata' contians the corresponding turbulent kinetic energy data
	# 'ukdata' holds the lower (FIRST COLUMN) and upper (SECOND COLUMN) bounds of k
	# 'veldata' contains the corresponding velocity component values
	# 'uveldata' holds the standard deviation of velocity in the x, y, and z
	##########
	posdata = np.zeros([nRows, 2])
	veldata = np.zeros([nRows, 3])
	uveldata = np.zeros([nRows, 3])
	kdata = np.zeros(nRows)
	ukdata = np.zeros([nRows, 2])
	print('  *Assuming 95% confidence means exactly 2 standard deviations (which is actually 95.45%)')
	# Extract the needed position and velocity data
	fmt = '%+.5e'
	for i in range(nRows):
		posdata[i,:] = np.array([matdata[i,1], matdata[i,2]])
		kdata[i] = np.array([matdata[i,18]])
		ukdata[i,:] = np.array([matdata[i,19], matdata[i,20]])
		veldata[i,:] = np.array([matdata[i,5], matdata[i,7], matdata[i,9]])
		uveldata[i,:] = np.array([matdata[i,6], matdata[i,8], matdata[i,10]])/2
		if (STATUS == True):
			if (i == 0):
				print("Data from {} includes:".format(projectPATH+'/data/'+data))
			print("Position (y,z): {}, {}".format(fmt%posdata[i,0],fmt%posdata[i,1]))
			print("\tVelocity x: {} +- {}".format(fmt%veldata[i,0],fmt%uveldata[i,0]))
			print("\tVelocity y: {} +- {}".format(fmt%veldata[i,1],fmt%uveldata[i,1]))
			print("\tVelocity z: {} +- {}".format(fmt%veldata[i,2],fmt%uveldata[i,2]))
			print("\tk (lower,upper): {} ({}, {})".format(fmt%kdata[i],fmt%ukdata[i,0],fmt%ukdata[i,1]))
			if (i == nRows-1):
				print("There are a total of {} inlet data points.".format(i))
	print("Done!")
	#####################################################################

	##########################################################################
	# Print the interpolated inlet conditions to the input file
	##########################################################################
	# Reformat the data points taken from OpenFOAM to work nicely with SciPy's 'griddata'
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
	# ### This integrates over the cross-section to determine the average velocity
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

	print("Gathering current U inlet conditions from '{}'... ".format(projectPATH+'/0/'+'U'), end='')
	lines = readin(projectPATH+'/0/'+'U')
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

	print("Writing interpolated U inlet conditions to '{}'... ".format(projectPATH+'/0/'+'U'), end='')
	# Print everything to a new file
	# Print the beginning of the file up till the inlet patch conditions
	fmt = '%.15f'
	f = open(projectPATH+'/0/'+'U', 'w')
	for i in np.arange(topStartPos+1):
		f.write(lines[i])
	# Write in the new inlet data
	f.write('        type            fixedValue;\n')
	f.write('        value           nonuniform List<vector>\n')
	f.write('{}\n'.format(nTopCells))
	f.write('(\n')
	for i in range(nTopCells):
		f.write('({} {} {})\n'.format(fmt%(newVx[i]), fmt%(newVy[i]), fmt%(newVz[i])))
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
		f.write('({} {} {})\n'.format(fmt%(newVx[i]), fmt%(newVy[i]), fmt%(newVz[i])))
	f.write(')\n')
	f.write(';\n')
	# Finish the rest of the file
	for i in np.arange(botEndPos,len(lines)):
		f.write(lines[i])
	f.close()
	print("Done!")
	##########################################################################

	##########################################################################
	# Repeat the exact same process above for the inlet turbulent kinetic energy values
	##########################################################################
	print("Interpolating k inlet conditions using '{}' scheme... ".format(scheme), end='')
	newK = griddata(posdata,
		kdata[:],
		matpos, method=scheme)
	print("Done!")

	print("Gathering current k inlet conditions from '{}'... ".format(projectPATH+'/0/'+'k'), end='')
	lines = readin(projectPATH+'/0/'+'k')
	topFound = False
	botFound = False
	for i, line in enumerate(lines):
		if (line.strip() == patches[0]):
			topStartPos = i+1
			topFound = True
		elif (topFound == True and line.strip() == '}'):
			topEndPos = i
			topFound = False
		elif (line.strip() == patches[1]):
			botStartPos = i+1
			botFound = True
		elif (botFound == True and line.strip() == '}'):
			botEndPos = i
			botFound = False
	print("Done!")

	print("Writing interpolated k inlet conditions to '{}'... ".format(projectPATH+'/0/'+'k'), end='')
	f = open(projectPATH+'/0/'+'k', 'w')
	for i in np.arange(topStartPos+1):
		f.write(lines[i])
	f.write('        type            fixedValue;\n')
	f.write('        value           nonuniform List<scalar>\n')
	f.write('{}\n'.format(nTopCells))
	f.write('(\n')
	for i in range(nTopCells):
		f.write('{}\n'.format(fmt%(newK[i])))
	f.write(')\n')
	f.write(';\n')
	for i in np.arange(topEndPos,botStartPos+1):
		f.write(lines[i])
	f.write('        type            fixedValue;\n')
	f.write('        value           nonuniform List<scalar>\n')
	f.write('{}\n'.format(nTopCells))
	f.write('(\n')
	for i in range(nTopCells,nCells):
		f.write('{}\n'.format(fmt%(newK[i])))
	f.write(')\n')
	f.write(';\n')
	for i in np.arange(botEndPos,len(lines)):
		f.write(lines[i])
	f.close()
	print("Done!")
	##########################################################################

	if (PLOT == True):
		plotData = np.array(['matvals', 'posdata', 'veldata', 'matpos',
			'kdata', 'newVx', 'newVy', 'newVz', 'newK'])
		for f in plotData:
			np.savetxt('inletPlotData/'+f, eval(f), delimiter=',')
		# Runs a python script in the background to plot the data
		os.system('python3 inletPlot.py &')