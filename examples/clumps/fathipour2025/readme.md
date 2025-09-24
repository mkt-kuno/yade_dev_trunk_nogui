# Algorithms for generating clumps (i.e. multi-sphere aggregates)

The adjacent `fathipour2025.py` file contains Python implementations of algorithms for generating clumps, or multi-sphere aggregates, from binary images of particle shape, as proposed by [Fathipour-Azar & Duriez (2025)](https://link.springer.com/article/10.1007/s11831-025-10256-1).

# Prerequisites

Usage requires to have Python 3 and the following libraries installed:

    csv
    trimesh
    numpy
    scipy
    scikit-learn

You can install these libraries using `pip` by running the following command:

    pip install csv trimesh numpy scipy scikit-learn

Note that installing `scikit-learn` provides access to a Python module named `sklearn` (as discussed, e.g., [therein](https://stackoverflow.com/q/38338822)) and that the corresponding installation, likewise to the cases of `numpy` and `scipy` is also possible through an `apt` package

```bash
sudo apt install python3-numpy python3-scipy python3-sklearn # to avoid numpy, scipy and scikit-learn in the above pip command
```

# Example

The following example illustrates possible commands to create a clump version of a superellipsoid particle (as initially described in a `LevelSet` fashion and a corresponding binary image):

```python
## Initial import considerations
import sys
sys.path.append('./')
import fathipour2025
from fathipour2025 import *
from yade import ymport
## For creating a binary image from a LevelSet-shaped particle
voxel_size = 0.1
lsB_id = O.bodies.append(levelSetBody('superellipsoid',extents = (1,2,0.5),epsilons = (0.4,1.2), spacing = voxel_size))
## We are now ready to actually work on a clump version of the above particle
csvFile = 'clump_members_se.csv'
## Following line is to use the Volume Coverage Algorithm, § 2.1.2 of the manuscript
volCov( numpy.array(O.bodies[lsB_id].shape.binarize()) # the binary image input
	,voxel_size # the underlying voxel size in the binary image
	,threshold = 0 # for keeping all clump members (obtained from the binary image voxels), whatever their radius
	,coverage = 1 # for covering all the particle volume, avoiding any a-priori volume simplification
	,out_csv = csvFile # in which file storing the clump members data
	,verbose = 0)
## Following line would be to use the K‐means Clustering Algorithm, § 2.1.3 of the manuscript and 250 clump members
# clustKMeans(numpy.array(O.bodies[lsB_id].shape.binarize()), voxel_size = voxel_size, n_clusters = 250, out_csv = csvFile)
## Actual YADE definition of the clump counterpart obtained in the above algorithm(s):
O.bodies.appendClumped(ymport.text(csvFile))
yade.qt.View() # where the observed position offset is because the Aabb of the clump version always starts at the space origin
```

# Prior development history

Development history prior to the Fall 2025 inclusion into YADE can be found at https://forge.inrae.fr/hadi.fathipour-azar/clumpGen

