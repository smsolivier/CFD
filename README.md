# CFD
NUEN 489 GEMIX Project

## run 
Python script to run interFoam in parallel. 
Does blockMesh and decomposePar before running and reconstructPar after running. 
Has command line options for 
* parallel decomposition (how many processors in x, y, z)
* number of volumes (split into x_inlet, x_mixing, y and z) 
* turbulence model (k-Epsilon or k-Omega) 
* switch for running in 3D 
* select the inlet velocity data (0.6 or 1)
* select the method used for inlet data interpolation (nearest|linear|cubic)
     *** 'linear' is reccommended (and default) for now, plan to add a b-spline interpolation method

See ./run -h for help. 

## Running on ada
* load Python/3.5.2 
* edit jobFile for wall time, number of processors, memory usage 
* run bsub < jobFile to submit job 

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

## clean 
bash script to remove OpenFOAM generated files. 

## Current Problems
* boundary conditions for omega, currently 1e7 at walls
* boundary conditions for k 
* inlet conditions from data 
* top and bottom inlets not mixing 
* flow is currently not turbulent 
