#! /bin/bash
cd ..
rm -rf Yade-OpenFOAM-coupling

if [ -z "$WM_PROJECT_VERSION" ]; then
    if [ -f "/root/OpenFOAM/OpenFOAM-v1906/etc/bashrc" ]; then
        bashrcPath=/root/OpenFOAM/OpenFOAM-v1906/etc/bashrc # compiled version (ubuntu18.04)
    elif [ -f "/usr/lib/openfoam/openfoam2312/etc/bashrc" ]; then
        bashrcPath=/usr/lib/openfoam/openfoam2312/etc/bashrc # precompiled package (ubuntu22.04)
    else #assume OFOAM6, use older coupling code
        bashrcPath=/root/OpenFOAM/OpenFOAM-6/etc/bashrc
    fi
    source $bashrcPath
fi


cp -rf ../../trunk/pkg/openfoam/coupling Yade-OpenFOAM-coupling
cd Yade-OpenFOAM-coupling
python3 setup.py

#### testing icoFoamYade ####
cd ../../../trunk/examples/openfoam/example_icoFoamYade
echo `pwd`
blockMesh
decomposePar
mkdir yadep
mkdir spheres
mpirun --allow-run-as-root -n 2 ../../../../install/bin/yade-ci scriptMPI.py

#### testing pimpleFoamYade ####
cd ../example_pimpleFoamYade
echo `pwd`
blockMesh
decomposePar
mkdir yadep
mkdir spheres
mpirun --allow-run-as-root -n 2 ../../../../install/bin/yade-ci scriptMPI.py
