#! /bin/bash
blockMesh
cp -r 0_org 0
decomposePar
mkdir spheres
# 6. Create a symbolic link to Yade Install
#ln -s /path/to/yade/install/bin/yade-exec yadeimport.py

mpirun --allow-run-as-root -n 2 python3 scriptMPI.py > log&
