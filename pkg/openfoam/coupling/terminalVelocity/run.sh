#! /bin/bash

touch yade.foam     # paraview reader

blockMesh
cp -r 0_org 0

setFields

decomposePar

mkdir spheres data
# 6. Create a symbolic link to Yade Install
#ln -s /path/to/yade/install/bin/yade-exec yadeimport.py

#In yade serial:
# python3 scriptMPI.py > log&
mpiexec -n 2 ../../../../../install/bin/yade-2024-09-15.git-382de81 scriptMPI.py
#In yade parallel
#mpirun --allow-run-as-root -n 2 python3 scriptMPI.py > log&
#
#
#
# restore0Dir
#
#
# runApplication blockMesh
#
# runApplication setFields
#
# runApplication decomposePar
#
# python3 boxWithInterfaceAndParticles_yade_script.py
# # mpirun --allow-run-as-root -n 2 python3 standingTank_script.py

#------------------------------------------------------------------------------
