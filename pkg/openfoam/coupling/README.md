# Yade-OpenFOAM-coupling
An OpenFOAM solver for realizing CFD-DEM simulations with the Open Source Discrete Element Solver Yade-DEM. 
 * Fast mesh search (including range based search) based on k-d Tree algorithm, faster than the original Octree search offered by OpenFOAM (mesh.findCell,  mesh.findNearest).
 * Gaussian interpolation of field variables. 
 * Simple point-force coupling (icoFoamYade) solver
 * Full 4-way coupling (pimpleFoamYade) solver (under validation).
 * Documentation : https://yade-dev.gitlab.io/trunk/FoamCoupling.html
 * Examples : https://gitlab.com/yade-dev/trunk/tree/master/examples/openfoam


## Build
* Source OpenFOAM's etc/bashrc
* Run setup.py : 
  * ``python3 setup.py``

## Running 
