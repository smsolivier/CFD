# CFD
NUEN 489 GEMIX Project

## run 
Python script to run interFoam in parallel with 4 cores. 
Does blockMesh and decomposePar before running and reconstructPar after running. 
Has command line options for 
	* parallel decomposition (how many processors in x, y, z)
	* number of volumes (split into x_inlet, x_mixing, y and z) 
	* turbulence model (k-Epsilon or k-Omega) 
	* switch for running in 3D 
Run as ./run Nx Ny Nz. See ./run -h for help. 

## kEpsilon/kOmega
Contains fvSchemes, fvSolution, turbulenceProperties files corresponding to k-Epsilon or k-Omega. 
run copies them into their correct place in system or constant. 

## mesh
Location of blockMeshDict. Top is the top half of the blockMeshDict that is common to both 3D and 2D runs. 2D and 3D change the boundary conditions to be appropriate. 
run alters local variables Nx1, Nx2, Ny, Nz which control the number of volumes used. 

## plotTimeStep
Plots the time step output to make sure the adjustable time step isn't going too small. 
Run as ./plotTimeStep or python3 plotTimeStep 
