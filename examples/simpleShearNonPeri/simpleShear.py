# J. Duriez jerome.duriez@inrae.fr
# Script to define a numerical model for a simple shear box (in a non-periodic space)
# In order to illustrate various Kinem...Engine which allow to perform different loadings on the box
# There is here one § for each Kinem...Engine, comment/uncomment them to observe what you want

#NB : the run of the script is paused at each plot window (so that there is time to observe it). Type "Return" in the Yade terminal to resume

from yade import plot
from yade.pack import *

#####################
# Material properties
#####################

O.materials.append(FrictMat(density=2600, young=4.0e9, poisson=.04, frictionAngle=.6))

########
# Bodies
########

#Definition of bodies in the numerical model: first, six boxes corresponding to sides of the simple shear box, then, the included particle sample

# Dimensions of the box:
length = 0.1
height = 0.02
width = 0.04
thickness = 0.001
# The boxes themselves:
leftBox = box(center=(-thickness / 2.0, (height) / 2.0, 0), extents=(thickness / 2.0, 5 * (height / 2.0 + thickness), width / 2.0), fixed=True, wire=True)
lowBox = box(center=(length / 2.0, -thickness / 2.0, 0), extents=(length / 2.0, thickness / 2.0, width / 2.0), fixed=True, wire=True)
rightBox = box(
        center=(length + thickness / 2.0, height / 2.0, 0), extents=(thickness / 2.0, 5 * (height / 2.0 + thickness), width / 2.0), fixed=True, wire=True
)
upBox = box(center=(length / 2.0, height + thickness / 2.0, 0), extents=(length / 2.0, thickness / 2.0, width / 2.0), fixed=True, wire=True)
behindBox = box(
        center=(length / 2.0, height / 2.0, -width / 2.0 - thickness / 2.0),
        extents=(2.5 * length / 2.0, height / 2.0 + thickness, thickness / 2.0),
        fixed=True,
        wire=True
)
inFrontBox = box(
        center=(length / 2.0, height / 2.0, width / 2.0 + thickness / 2.0),
        extents=(2.5 * length / 2.0, height / 2.0 + thickness, thickness / 2.0),
        fixed=True,
        wire=True
)
O.bodies.append([leftBox, lowBox, rightBox, upBox, behindBox, inFrontBox])

# The particle packing
expected_poros = 0.76  # expected value of porosity after calling makeCloud below
nSpheres = 2000
rMean = pow(((1 - expected_poros) * length * height * width) / (nSpheres * 4.0 / 3.0 * pi), 1.0 / 3.0)
sp = SpherePack()
sp.makeCloud(Vector3(0, 0.0, -width / 2.0), Vector3(length, height, width / 2.0), rMean, .15)
sp.toSimulation()
# print('Obtained porosity is',porosity(),'while mean radius had been computed on the assumption of getting',expected_poros)

#########
# Engines
#########

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb()]),
        InteractionLoop([Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()], [Ip2_FrictMat_FrictMat_FrictPhys()], [Law2_ScGeom_FrictPhys_CundallStrack()]),
        NewtonIntegrator(damping=.1),
        PyRunner(iterPeriod=50, command='defData()'),
        KinemCTDEngine(compSpeed=0.4, targetSigma=50e3)  # for defining the initial compression
        ,
        PyRunner(iterPeriod=1, command='prepareSwitchToShear(nCycShear)', label='checker')
]
nCycShear = 20000  # 20000 will be the number of active shear iterations, when time will come


def defData():
	# NB: mechanical behavior is interpreted in an interface manner with displacement variables
	plot.addData(
	        fy=O.forces.f(3)[1],  # vertical component of the force sustained by the upper side of the shear box
	        fx=O.forces.f(3)[0],  # horizontal component of the force sustained by the upper side of the shear box
	        step=O.iter,
	        gamma=O.bodies[3].state.pos[0] - length / 2.0,  # relative shear displacement
	        u=-O.bodies[3].state.pos[1] + (height + thickness / 2.0)  # relative normal displacement (positive in compression)
	)


def prepareSwitchToShear(nCycShear):
	if O.iter > 82000:
		O.pause()
		O.engines = O.engines[:5] + [KinemCNDEngine(shearSpeed=(length / 7.0) / (nCycShear * O.dt), gammalim=length / 7.0)] + [O.engines[-1]]
		checker.command = 'prepareEndOfLoading(nCycShear)'


def prepareEndOfLoading(nCycShear):
	if O.iter - 82000 > 1.15 * nCycShear:
		print('Stopping')
		O.pause()


O.dt = .4 * PWaveTimeStep()

###########################
# Simulation will now run ! With various loading phases
###########################

from yade import qt

qt.View()

#---- [Oedometer-like] Compression with KinemCTDEngine defined in the above O.engines list ----

print('\nRunning in progress: the shear box is being compressed')
O.run()
plot.plots = {'step': ('fy'), 'step ': ('u')}
plot.plot()
print('Plotting curves (with fy = normal force and u = normal displacement). Type Return to go ahead once things have stopped to move (after 82000 steps)\n')
input()

#---- Shear at constant normal displacement: KinemCNDEngine (defined in prepareSwitchToShear()) ----

print('The shear box is now being sheared')
O.run()
plot.plots = {'step': ('gamma'), 'gamma': ('fx', 'fy')}
plot.plot()
print(
        'Plotting curves (gamma = tangential displacement,  fx = tangential force)) during shear loading. Will stop by itself after a bit more than', nCycShear,
        'additional iterations\n'
)

# #---- A re-compression, from this initial sheared state: KinemCTDEngine again ----
# O.engines=O.engines[:5]+[KinemCTDEngine(compSpeed=0.5,sigma_save=(),temoin_save=(),targetSigma=80000.0,LOG=False)]
# print('Be patient, running in progress (the sample is being again compressed, from this sheared state)')
# O.run(10000,True)
# plot.plots={'u':('fx','fy',)}
# plot.plot(subPlots=False)
# print('Plotting curve. Type Return to go ahead\n')
# input()

##---- Shear at constant normal load/stress ----
#nCycShear = 20000
#O.engines=O.engines[:6]+[KinemCNLEngine(shearSpeed=(length/10.0)/(nCycShear*O.dt),gamma_save=(),temoin_save=(),gammalim=length/10.0,LOG=False)]
#O.run(int(1.15*nCycShear),True)
#plot.plots={'step':('gamma','u',)}
#plot.plot()
#input()
#plot.plots={'gamma':('fx','fy',)}
#plot.plot()
#input()

#---- Shear at constant normal stifness ----
#nCycShear = 20000
#O.engines=O.engines[:6]+[KinemCNSEngine(shearSpeed=(length/10.0)/(nCycShear*O.dt),gammalim=length/10.0,LOG=False,KnC=1)]
#O.run(int(1.15*nCycShear),True)
#plot.plots={'step':('gamma','u',)}
#plot.plot()
#input()
#plot.plots={'gamma':('fx','fy',)}
#plot.plot()
#input()
#plot.plots={'u':('fy',)}
#plot.plot()
#input()
