#! /bin/bash

cd ..
rm -rf Yade-OpenFOAM-coupling

if [ -z "$WM_PROJECT_VERSION" ]; then
    if [ -f "/root/OpenFOAM/OpenFOAM-v1906/etc/bashrc" ]; then
        bashrcPath=/root/OpenFOAM/OpenFOAM-v1906/etc/bashrc # compiled version (ubuntu18.04)
        git clone -b supportOFoam1906 https://gitlab.com/yade-dev/Yade-OpenFOAM-coupling.git
    elif [ -f "/usr/lib/openfoam/openfoam2312/etc/bashrc" ]; then
        bashrcPath=/usr/lib/openfoam/openfoam2312/etc/bashrc # precompiled package (ubuntu22.04)
        git clone -b supportOFoam1906 https://gitlab.com/yade-dev/Yade-OpenFOAM-coupling.git
    else #assume OFOAM6, use older coupling code
        bashrcPath=/root/OpenFOAM/OpenFOAM-6/etc/bashrc
        git clone -b supportOFoam1906 https://gitlab.com/yade-dev/Yade-OpenFOAM-coupling.git
    fi
    source $bashrcPath
else  # some bashrc have been sourced already (e.g. in .gitlab-ci.yml for the pipeline)
    git clone -b supportOFoam1906 https://gitlab.com/yade-dev/Yade-OpenFOAM-coupling.git
fi

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

mpirun --allow-run-as-root -n 4 ../../../install/bin/yade-ci scriptMPI.py

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

mpirun --allow-run-as-root -n 4 ../../../install/bin/yade-ci scriptMPI.py
