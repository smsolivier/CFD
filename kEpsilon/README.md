# CFD
NUEN 489 GEMIX Project

## run 
Python script to run interFoam in parallel with 4 cores. 
Does blockMesh and decomposePar before running and reconstructPar after running. 
Run as ./run Nx Ny Nz. See ./run -h for help. 

## plotTimeStep
Plots the time step output to make sure the adjustable time step isn't going too small. 
Run as ./plotTimeStep or python3 plotTimeStep 
