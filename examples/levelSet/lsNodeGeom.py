# A y-axis bouncing (and x-axis translating) LS-sphere on a LS-floor to illustrate the use of LSnodeGeom

##############
# The bodies #
##############

# The LS floor (with manual definition and assignment of distance field):
xMax, yMax, zMax, spacing = 5, 1, 2, 1
nGP = (numpy.array([xMax, yMax, zMax]) * 2 + 1) / spacing
nGP = Vector3i(int(nGP[0]), int(nGP[1]), int(nGP[2]))
grid = RegularGrid(min=Vector3(-xMax, -yMax, -zMax), spacing=spacing, nGP=nGP)
X, Y, Z = numpy.meshgrid(numpy.linspace(-xMax, xMax, nGP[0]), numpy.linspace(-yMax, yMax, nGP[1]), numpy.linspace(-zMax, zMax, nGP[2]), indexing='ij')
distField = Y
O.bodies.append(
        levelSetBody(grid=grid, distField=distField.tolist(), dynamic=False)
)  # because we know here that surface nodes will not express on that LS floor (bigger than the sphere), we do not need to pay attention to those and n_neighborsNodes

# The LS sphere:
lsSph = O.bodies.append(
        levelSetBody('sphere', (-0.9 * xMax, 1.9, 0), 0.5, spacing=0.05, nSurfNodes=1602, n_neighborsNodes=6)
)  # mind the n_neighborsNodes assignment
lsSph = O.bodies[lsSph]
sphMass = lsSph.state.mass
lsSph.state.vel = Vector3(1, 0, 0)

###############
# The engines #
###############

grav = 10
kn = sphMass * grav / 0.001
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_LevelSet_Aabb()],
                              verletDist=0)  # absence of Sphere bodies prevents an automatic determination of verletDist. We use here 0 (non-optimal)
        ,
        InteractionLoop(
                [Ig2_LevelSet_LevelSet_LSnodeGeom()],
                [Ip2_FrictMat_FrictMat_FrictPhys(kn=MatchMaker(algo='val', val=kn), ks=MatchMaker(algo='val', val=kn))],
                [Law2_ScGeom_FrictPhys_CundallStrack(sphericalBodies=False, neverErase=True)]  # mind the neverErase = True
                ,
                warnRoleIg2=False
        )  # mind the warnRoleIg2 = False
        ,
        NewtonIntegrator(gravity=Vector3(0, -grav, 0), damping=0.05),
        VTKRecorder(recorders=["lsBodies"], iterPeriod=50, multiblockLS=True),
        PyRunner(iterPeriod=1, command='saveData()')
]
from yade import plot


def saveData():
	if O.interactions.has(0, 1, True):
		cont = O.interactions[0, 1]
		nodeIdx = cont.geom.surfNodeIdx
		un = cont.geom.penetrationDepth
	else:
		nodeIdx = nan
		un = nan
	plot.addData(it=O.iter, pos=lsSph.state.pos, velLin=lsSph.state.vel, velAng=lsSph.state.angVel, un=un, nodeIdx=nodeIdx)


O.dt = 0.5 * (sphMass / kn)**0.5
plot.plots = {'it': ('pos_x', 'pos_y', 'pos_z'), ' it': ('velAng_x', 'velAng_y', 'velAng_z'), 'it ': 'un', ' it ': 'nodeIdx'}
plot.plot()
O.saveTmp()
yade.qt.View()
nIter = 3000
O.run(nIter)
# post-processing in Paraview with, e.g., pvVisu(itSnap = [50*it for it in range(int(nIter/50))], multiblockLS = True)
