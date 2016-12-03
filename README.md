# CFD
NUEN 489 GEMIX Project
## To Do 
* Update inlet to write boundary conditions for omega, epsilon using formulas 
* simpleFoam comparison to make sure interFoam isn't buggy 
* Sensitivity analysis
* Model error 

## run 
Python script to run interFoam in parallel. 
Does blockMesh and decomposePar before running and reconstructPar after running. 
Has command line options for 
* parallel decomposition (how many processors in x, y, z)
* number of volumes (split into x_inlet, x_mixing, y and z) 
* turbulence model (k-Epsilon or k-Omega) 
* switch for running in 3D 
* select which case to run 
* select the method used for inlet data interpolation (nearest|linear|cubic)
    'linear' is recommended (and default) for now, plan to add a b-spline interpolation method

See ./run -h for help. 

## Running on ada
* edit jobFile for wall time, number of processors, memory usage 
* edit ./run command to include proper flags (don't plot anything on compute nodes) 
* run bsub < jobFile to submit job 
* bjobs to make sure submitted 
* bpeek to check stdout and stderr 

## kEpsilon/kOmega
Contains fvSchemes, fvSolution, turbulenceProperties files corresponding to k-Epsilon or k-Omega. 
./run copies them into their correct place in system or constant before running interFoam. 

## mesh
Location of blockMeshDict. Top is the top half of the blockMeshDict that is common to both 3D and 2D runs. 2D and 3D change the boundary conditions to be appropriate for 2 and 3 D geometry. simple is the blockMeshDict for 3D geometry without angled inlets. 

./run alters local variables Nx1, Nx2, Ny, Nz which control the number of volumes used. 

## plotTime.py
Plots the time step output to make sure the adjustable time step isn't going too small. 
Run as ./plotTime.py or python3 plotTime

## generateInlet.py 
Reads inlet data and writes to 0/U and 0/k before running. 

## inletPlot.py
Used along with the -PLOT flag during the program execution. Shows a 3D scatter plot of the cell center locations overlayed on the upper and lower inlets and a color map of the 4 inlet conditions (x, y, and z velocity and turbulent kinetic energy) for each cell. This allows to visually verify that the cell centers are in the correct position and appear to take the correct distribution.

## fRe.py
Reads in a CSV file of the axial pressure to calculate the fRe product. 

## readFoam.py 
Reads latest time directory and writeCellCentres output to plot the velocity profile against the experimental data. 

./readFoam.py case_number 

## readParaView.py 
Utility to read output csv's. 

Example use:

import readParaView as rpv 

df, names = rpv.read('output.csv', 'U:0', 'U:1', 'U:2', 'Points:1')

returns the columns of data corresponding to U:0, U:1, U:2, Points:1 and a list of the header names. 

## clean 
bash script to remove OpenFOAM generated files.  
