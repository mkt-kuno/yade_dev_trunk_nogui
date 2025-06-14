# Script to simulate a particle in a box containing two fluids
# 2024 (c) by Vasileios Angelidakis <v.angelidakis@qub.ac.uk>
# Based on the script CollapseExample.py of Eduard Puig Montella

import os, shutil
from yadeimport import *
from yade import mpy as mp
from yade import ymport, geom

parallelYade = False  #mpirun --allow-run-as-root -n 2 python3 scriptMPI.py , if False  python3 scriptMPI.py
numProcOF = 2

NSTEPS = 2000000
saveVTK = 1
radius = 0.005
# --------------------------------------------------------------------------------------------------
# Create materials
#O.materials.append(ViscElMat(en=0.9, et=0.9, young=2.11e7, poisson=0.3, density=7980, frictionAngle=atan(0.52), label='spheres'))
#O.materials.append(ViscElMat(en=0.9, et=0.9, young=2.11e7, poisson=0.3, density=0   , frictionAngle=atan(0.52), label='walls'    ))
O.materials.append(FrictMat(young=2.1e7, poisson=0.30, density=200, frictionAngle=atan(0.52), label='spheres'))
O.materials.append(FrictMat(young=1.0e7, poisson=0.23, density=7980, frictionAngle=atan(0.43), label='walls'))

# --------------------------------------------------------------------------------------------------
# Generate spherical particles

O.bodies.append(sphere((0.5, 0.5, 0.2), radius, material='spheres'))
O.bodies[0].state.vel = [0, 0, 0]

# sphereIDs = [b.id for b in O.bodies if type(b.shape) == Sphere] # FIXME: Amend this to count clumps
print("Number of particles:", len(O.bodies))

# --------------------------------------------------------------------------------------------------
# Set up OpenFOAM coupling parameters
fluidCoupling = FoamCoupling()
fluidCoupling.couplingModeParallel = parallelYade
fluidCoupling.isGaussianInterp = True
#use pimpleFoamYade for gaussianInterp (only in serial mode)
sphereIDs = [b.id for b in O.bodies if type(b.shape) == Sphere]
#wallsLockId = O.bodies.append(box(center= (L_tank, H/2.,W/2.),extents=(0,H/2.,W/2.),fixed=True,wire=False,color = (1.,0.,0.),material='walls'))
'''The yade specific (icoFoamYade, pimpleFoamYade) OpenFOAM solver can be found in $FOAM_USER_APPBIN, (
# full path here, the scond argument, 2 is the number of FoamProcs. '''
# fluidCoupling.SetOpenFoamSolver(os.environ.get('FOAM_USER_APPBIN')+'/icoFoamYade', 2)
# it also work without path after sourcing OFoam's bashrc
fluidCoupling.SetOpenFoamSolver("interFoamYadev2312", numProcOF)

# --------------------------------------------------------------------------------------------------
# Define engines
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()], label='collider'),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
                [Ip2_FrictMat_FrictMat_MindlinPhys(en=0.9, es=0.9)],  #Ip2_ViscElMat_ViscElMat_ViscElPhys
                [Law2_ScGeom_MindlinPhys_Mindlin()],  #Law2_ScGeom_ViscElPhys_Basic
                label="InteractionLoop"
        ),
        #	GlobalStiffnessTimeStepper(timestepSafetyCoefficient=0.7, timeStepUpdateInterval=100, parallelMode=parallelYade, label="ts"),
        fluidCoupling,  #to be called after timestepper
        NewtonIntegrator(gravity=(0, 0, -9.81), damping=0.0, label='newton'),  #(-0,-9.81,0)
        PyRunner(command='printPosition()', iterPeriod=saveVTK),
        VTKRecorder(fileName='results/spheres/3d-vtk-', recorders=['spheres'], parallelMode=parallelYade, iterPeriod=saveVTK)
]

# Example usage
#O.dt=0.20*RayleighWaveTimeStep(); print(O.dt*10e6) # FIXME: Check this
# O.dt=2.13e-6; # FIXME: Check how this compares with 0.2 of the Rayleigh critical timestep
O.dt = 1e-6


def printPosition():
	print(O.bodies[0].state.pos)


#collider.verletDist = 0.0075
mp.YADE_TIMING = False
mp.FLUID_COUPLING = True
mp.VERBOSE_OUTPUT = False
mp.USE_CPP_INTERS = True
mp.ERASE_REMOTE_MASTER = True
mp.REALLOC_FREQUENCY = 0
mp.fluidBodies = sphereIDs
#mp.commSplit = True
mp.DOMAIN_DECOMPOSITION = True

########### sedimentation ###########
# #
#while 1:
#	mp.mpirun(2000)
#	unb = unbalancedForce()
#	if unb < 0.01 :
#		newton.damping=0.1
#		for b in O.bodies:
#			if b.id ==wallsLockId:
#				mp.bodyErase(b.id)
#				print ("wallsLockId removed!")
#		break

#dataFProfile.dead=0

mp.mpirun(NSTEPS)
mp.mprint("RUN FINISH")
# fluidCoupling.killMPI()
exit()

## --------------------------------------------------------------------------------------------------
## Visualise the scene
#from yade import qt
#v=qt.View()
#v.ortho=True

#v.eyePosition=Vector3(0,-2*H,H/2)
#v.upVector=Vector3(0,0,1)
#v.viewDir=Vector3(0,1,0)

#rndr=qt.Renderer()

#O.saveTmp()
## --------------------------------------------------------------------------------------------------
## Run the script serially
#O.run(10000)
