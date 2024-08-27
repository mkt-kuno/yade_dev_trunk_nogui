#! /bin/bash
blockMesh
cp -r 0_org 0
decomposePar
mpirun --allow-run-as-root -n 2 ../../../../install/bin/yade-ci scriptMPI.py
