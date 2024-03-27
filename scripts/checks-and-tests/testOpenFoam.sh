#! /bin/bash

cd ..
rm -rf Yade-OpenFOAM-coupling


if [ -f "/root/OpenFOAM/OpenFOAM-v1906/etc/bashrc" ]; then
    bashrcPath=/root/OpenFOAM/OpenFOAM-v1906/etc/bashrc # compiled version (ubuntu18.04)
    git clone -b supportOFoam1906 https://gitlab.com/yade-dev/Yade-OpenFOAM-coupling.git
elif [ -f "/usr/lib/openfoam/openfoam2012/etc/bashrc" ]; then
    bashrcPath=/usr/lib/openfoam/openfoam2012/etc/bashrc # precompiled package (ubuntu22.04)
    git clone -b supportOFoam1906 https://gitlab.com/yade-dev/Yade-OpenFOAM-coupling.git
else #assume OFOAM6, use older coupling code
    bashrcPath=/root/OpenFOAM/OpenFOAM-6/etc/bashrc
    git clone https://github.com/dpkn31/Yade-OpenFOAM-coupling.git
fi

source $bashrcPath
cd Yade-OpenFOAM-coupling
./Allclean
./Allwmake

cd ../trunk/examples/openfoam/example_icoFoamYade
echo `pwd`
blockMesh
decomposePar
mkdir yadep
mkdir spheres



if [ -f icoFoamYade ]; then
    echo 'File exists.'
else
    echo 'File does not exist.'
fi

mpirun --allow-run-as-root -n 4 ../../../install/bin/yade-ci scriptMPI.py


