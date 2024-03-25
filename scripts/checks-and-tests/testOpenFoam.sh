#! /bin/bash

ls -la /root/OpenFOAM/OpenFOAM-v1906/etc/bashrc
source  /root/OpenFOAM/OpenFOAM-v1906/etc/bashrc

cd ..
rm -rf Yade-OpenFOAM-coupling
if [ -v OFOAM6 ]; then
	git clone https://github.com/dpkn31/Yade-OpenFOAM-coupling.git
else
	git clone -b supportOFoam1906 https://gitlab.com/yade-dev/Yade-OpenFOAM-coupling.git
fi
cd Yade-OpenFOAM-coupling
# git checkout yadeTestPar #why?
./Allclean
./Allwmake

cd ../trunk/examples/openfoam/example_icoFoamYade
echo `pwd`
echo "---------- creating symlink yadeimport.py ----------"
ln -s ../../../install/bin/yade-ci ./yadeimport.py
ls -la ./yadeimport.py
blockMesh
decomposePar
mkdir yadep
mkdir spheres

if [ -f /root/OpenFOAM/-v1906/platforms/linux64GccDPInt32Opt/bin/icoFoamYade ]; then
    echo 'File exists.'
else
    echo 'File does not exist.'
fi


# mpiexec --allow-run-as-root -n 1 python3 scriptYade.py : -n 2 icoFoamYade -parallel
# mpiexec --allow-run-as-root -n 4 ../../../install/bin/yade-ci scriptMPI.py : -n 2 icoFoamYade -parallel
echo "_____ RUNING mpirun --allow-run-as-root -n 1 python3 scriptYade.py"
sleep 2
mpirun --allow-run-as-root -n 1 python3 scriptYade.py

echo "_____ RUNING mpirun --allow-run-as-root -n 4 ../../../install/bin/yade-ci scriptMPI.py"
sleep 2
mpirun --allow-run-as-root -n 4 ../../../install/bin/yade-ci scriptMPI.py
# echo "_____ RUNING ../../../install/bin/yade-ci scriptMPI.py"
# sleep 2
# ../../../install/bin/yade-ci scriptMPI.py

