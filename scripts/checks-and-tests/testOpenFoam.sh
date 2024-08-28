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

#git clone https://github.com/dpkn31/Yade-OpenFOAM-coupling.git
cp -rf trunk/pkg/openfoam/coupling Yade-OpenFOAM-coupling
cd Yade-OpenFOAM-coupling
python3 setup.py
#./Allclean
#./Allwmake

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

mpirun --allow-run-as-root -n 2 ../../../../install/bin/yade-ci scriptMPI.py

cd ../example_pimpleFoamYade
echo `pwd`
blockMesh
decomposePar
mkdir yadep
mkdir spheres

if [ -f pimpleFoamYade ]; then
   echo 'File exists.'
else
   echo 'File does not exist.'
fi

mpirun --allow-run-as-root -n 2 ../../../../install/bin/yade-ci scriptMPI.py
