# Hopf-Cole / Varadhan distance 
Compute an approximation of the distance function to the boundary of a domain by solving a screened Poisson equation and using Varadhan's formula / the Hopf-Cole transform. See equations (4) and (5) in this [paper](https://arxiv.org/abs/1204.6216) or equations (6) - (11) in this [one](https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.12611). 

## Examples 
Clone this repo, go to the directory 'src' and run one of the scripts 'demo_fig' or 'demo_to_vtk' within Matlab. 
The first script will produce figures using Matlab's builtins functions, while the second will generate VTK files (legacy, text format) that can be visualized with [paraview](https://www.paraview.org/). 

## Input data 
The scripts above expect as input a triangulated 2D domain. The file format exported by the program [triangle](https://www.cs.cmu.edu/~quake/triangle.html) is used. Examples of data for two simple domains (a square and a disk), and a more complicated one (a rider) are provided in the directory 'data'. 
