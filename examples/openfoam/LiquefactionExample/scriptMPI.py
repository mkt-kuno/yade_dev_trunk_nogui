

# -*- encoding=utf-8 -*-
import os
from yadeimport import *
from yade import mpy as mp
import numpy as np

import mpi4py.MPI as MPI
# path = "spheres"
# isExist = os.path.exists(path)
# if not isExist:
#    os.makedirs(path)
#    print("The new directory"+path+" is created!")

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

compFricCoef = 0.47 # initial contact friction during the confining phase
finalFricCoef = 0.47 # contact friction during the deviatoric loading
wallFricCoef = 0.19 # contact friction during the deviatoric loading
densitySpheres = 2500
g = 9.81

young=6.3e7
poissonR=0.2
L = 0.4
W = 0.04
H = 0.1
Radius=70e-6/2.*100
xStep=Radius*0.4
porosity=0.6
restitCoef=0.91

density = 1000
NSTEPS = 9000000
saveVTK=1000



# O.materials.append(ViscElMat(en=restitCoef, et=1., young=young, poisson=poissonR, density=densitySpheres, frictionAngle=compFricCoef, label='spheremat'))
# O.materials.append(ViscElMat(en=restitCoef, et=1., young=young, poisson=poissonR, density=0, frictionAngle=wallFricCoef, label='walls'))

O.materials.append(FrictMat(young=young,poisson=poissonR,frictionAngle=radians(compFricCoef),density=densitySpheres,label='spheres'))
O.materials.append(FrictMat(young=young,poisson=poissonR,frictionAngle=radians(25),density=0,label='walls'))

#
# mx = Vector3(L,H,W)
# mn = Vector3(-0.0, -0.0, -0.0)
# mnSph = Vector3(L, H, W)

mx = Vector3(L,H,W)
mn = Vector3(-0.0, -0.0, -0.0)
mnSph = Vector3(L, H, W)

sp = pack.SpherePack()
sp.makeCloud(mn, mnSph, rMean=Radius, rRelFuzz=0.0,porosity=porosity,periodic=False)  #"seed" make the "random" generation always the same
sp.toSimulation(material='spheres')
print ("Numer of particles:", len(O.bodies))



walls = aabbWalls([mn, mx], thickness=0,material='walls')
wallIds = O.bodies.append(walls)

fluidCoupling = FoamCoupling()
fluidCoupling.couplingModeParallel = True
fluidCoupling.isGaussianInterp = True
#use pimpleFoamYade for gaussianInterp (only in serial mode)
sphereIDs = [b.id for b in O.bodies if type(b.shape) == Sphere]
# for b in O.bodies:
# 	if type(b.shape) == Sphere:
# 		b.state.blockedDOFs='xyzXYZ'

'''The yade specific (icoFoamYade, pimpleFoamYade) OpenFOAM solver can be found in $FOAM_USER_APPBIN, (
# full path here, the scond argument, 2 is the number of FoamProcs. '''
# fluidCoupling.SetOpenFoamSolver(os.environ.get('FOAM_USER_APPBIN')+'/icoFoamYade', 2)
# it also work without path after sourcing OFoam's bashrc
fluidCoupling.SetOpenFoamSolver("pimpleFoamYade", 2)





newton=NewtonIntegrator(gravity=(-0,-g,0),damping=0.1)

O.engines=[
 ForceResetter(),
 InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb()], label='collider', allowBiggerThanPeriod=True),
 #InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb(),Bo1_Wall_Aabb()]),
	InteractionLoop(
	        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom(), Ig2_Wall_Sphere_ScGeom()], [Ip2_FrictMat_FrictMat_FrictPhys()],
	        [Law2_ScGeom_FrictPhys_CundallStrack()]
	),
        GlobalStiffnessTimeStepper(timestepSafetyCoefficient=0.7, timeStepUpdateInterval=100, parallelMode=True, label="ts"),
        fluidCoupling,  #to be called after timestepper
        newton,
        VTKRecorder(fileName='spheres/3d-vtk-', recorders=['spheres','boxes'], parallelMode=True, iterPeriod=saveVTK),
        PyRunner(command="UpPlot()", iterPeriod=saveVTK, label="UpPlotLabel"),
        PyRunner(iterPeriod=saveVTK,command='saveProfile(O.iter)',dead=1,label="dataFProfile")
]





def UpPlot():
	if rank==1:
		points=[]
		for b in O.bodies:
				if type(b.shape) == Sphere:
						points.append(b.state.pos[1])

		height = sum(points) / len(points)
		phi=1-voxelPorosity(200,Vector3(Radius,Radius,Radius),Vector3(L-Radius,6*Radius,W-Radius))
		with open('TimeData.txt','a') as f:
				f.write(str(O.iter) + " "+str(height) + " " +str(phi) + " "+"\n")
				f.close




collider.verletDist = 0.0075
mp.YADE_TIMING = False
mp.FLUID_COUPLING = True
mp.VERBOSE_OUTPUT = False
mp.USE_CPP_INTERS = True
mp.ERASE_REMOTE_MASTER = True
mp.REALLOC_FREQUENCY = 0
mp.fluidBodies = sphereIDs
#mp.commSplit = True
mp.DOMAIN_DECOMPOSITION = True
#
########## sedimentation ###########
#
# while 1:
# 	mp.mpirun(1000)
# 	unb = unbalancedForce()
# 	mp.mprint ("forces:", unb)
# 	with open('WallDeleted.txt', 'a') as f:
# 		 f.write("forces: " + str(unb)+ "\n")
# 		 f.close
# 	# if unb < 0.01 or O.iter >2e4:
# 	if O.iter >2e4:
# 		mp.mprint ("gravitational deposition completed")
# 		setContactFriction(finalFricCoef)
# 		newton.damping=0.05
#
# 		# with open('WallDeleted.txt', 'a') as f:
# 		#     f.write("Gate was removed at iteration: " + str(O.iter) + " and time: " +  str(O.time) + "\n")
# 		#     f.close
# 		# break

mp.mpirun(NSTEPS)
mp.mprint("RUN FINISH")
#fluidCoupling.killMPI()
exit()
