import os
from yadeimport import *
from yade import mpy as mp

# path = "spheres"
# isExist = os.path.exists(path)
# if not isExist:
#    os.makedirs(path)
#    print("The new directory"+path+" is created!")


compFricCoef = 0.0 # initial contact friction during the confining phase
finalFricCoef = 0.47 # contact friction during the deviatoric loading
wallFricCoef = 0.19 # contact friction during the deviatoric loading
densitySpheres = 2178
g = 9.81

parallelYade=False #mpirun --allow-run-as-root -n 2 python3 scriptMPI.py , if False  python3 scriptMPI.py
numProcOF=2

young=6.3e7
poissonR=0.2
viscosityF=0.0291
L_tank = 0.04
L = 0.4
W = 0.04
H = 0.325
Radius=3.84e-3/2.*2
porosity=0.6
restitCoef=0.91

density = 1000
NSTEPS = 2000000
saveVTK=1000

O.materials.append(ViscElMat(en=restitCoef, et=1., young=young, poisson=poissonR, density=densitySpheres, frictionAngle=compFricCoef, label='spheremat'))
O.materials.append(ViscElMat(en=restitCoef, et=1., young=young, poisson=poissonR, density=0, frictionAngle=wallFricCoef, label='walls'))


mx = Vector3(L,H,W)
mn = Vector3(-0.0, -0.0, -0.0)
mnSph = Vector3(L_tank, H*0.5, W)
# sp = pack.SpherePack()
# sp.makeCloud(mn, mnSph, rMean=Radius, rRelFuzz=0.0781,porosity=porosity,periodic=False)  #"seed" make the "random" generation always the same
sp=pack.randomDensePack(pack.inAlignedBox(mn,mnSph),radius=Radius,rRelFuzz=0.0781,returnSpherePack=True,seed=1)


sp.toSimulation(material='spheremat')
# sphereIDs = [b.id for b in O.bodies if type(b.shape) == Sphere]
print ("Numer of particles:", len(O.bodies))



walls = aabbWalls([mn, mx], thickness=0,material='walls')
wallIds = O.bodies.append(walls)

fluidCoupling = FoamCoupling()
fluidCoupling.couplingModeParallel = parallelYade
fluidCoupling.isGaussianInterp = True
#use pimpleFoamYade for gaussianInterp (only in serial mode)
sphereIDs = [b.id for b in O.bodies if type(b.shape) == Sphere]
wallsLockId = O.bodies.append(box(center= (L_tank, H/2.,W/2.),extents=(0,H/2.,W/2.),fixed=True,wire=False,color = (1.,0.,0.),material='walls'))


'''The yade specific (icoFoamYade, pimpleFoamYade) OpenFOAM solver can be found in $FOAM_USER_APPBIN, (
# full path here, the scond argument, 2 is the number of FoamProcs. '''
# fluidCoupling.SetOpenFoamSolver(os.environ.get('FOAM_USER_APPBIN')+'/icoFoamYade', 2)
# it also work without path after sourcing OFoam's bashrc
fluidCoupling.SetOpenFoamSolver("pimpleFoamYade", numProcOF)


def saveProfile(k):
    for b in O.bodies:
        if type(b.shape) == Sphere:
                with open('data/profile_%i.txt'%k, 'a') as f:
                            f.write(str(b.state.pos[0]) + " " +str(b.state.pos[1]) + " " + str(b.state.pos[2])+ "\n")
                            f.close
    return

newton=NewtonIntegrator(gravity=(-0,-g,0),damping=0.4)

O.engines=[
 ForceResetter(),
 InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb()], label='collider', allowBiggerThanPeriod=True),
 #InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb(),Bo1_Wall_Aabb()]),
 InteractionLoop(
  [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom(), Ig2_Wall_Sphere_ScGeom()],
  #[Ig2_Sphere_Sphere_ScGeom6D(interactionDetectionFactor=enlFactor),Ig2_Box_Sphere_ScGeom6D(interactionDetectionFactor=enlFactor), Ig2_Wall_Sphere_ScGeom()],
  [Ip2_ViscElMat_ViscElMat_ViscElPhys()],
  [Law2_ScGeom_ViscElPhys_Basic()],label="InteractionLoop"
 ),
        GlobalStiffnessTimeStepper(timestepSafetyCoefficient=0.7, timeStepUpdateInterval=100, parallelMode=parallelYade, label="ts"),
        fluidCoupling,  #to be called after timestepper
        newton,
        VTKRecorder(fileName='spheres/3d-vtk-', recorders=['spheres','boxes'], parallelMode=parallelYade, iterPeriod=saveVTK),
        PyRunner(iterPeriod=saveVTK,command='saveProfile(O.iter)',dead=1,label="dataFProfile")
]


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

########## sedimentation ###########
# #
while 1:
	mp.mpirun(2000)
	unb = unbalancedForce()
	if unb < 0.01 :
		newton.damping=0.1
		for b in O.bodies:
			if b.id ==wallsLockId:
				mp.bodyErase(b.id)
				print ("wallsLockId removed!")
		break


dataFProfile.dead=0

mp.mpirun(NSTEPS)
mp.mprint("RUN FINISH")
#fluidCoupling.killMPI()
exit()
