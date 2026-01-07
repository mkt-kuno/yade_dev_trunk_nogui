### constants
nTheta = 4
nPhi = 8
interParticleExchangeCoeff = 25

### material
matSeg1 = CohFrictMatSeg()

### bodies
s1_seg = sphere((0,0,0), 1, material = matSeg1, fixed = True )
s2_seg = sphere((2.001,0,0), 1, material = matSeg1, fixed = True )
s1_id_seg, s2_id_seg = O.bodies.append([s1_seg, s2_seg])

s2_seg.state.vel = (-1.0,0,0)

### initialize updater
updater = SegmentedStateUpdater(thetaResolution = nTheta, phiResolution = nPhi, interParticleExchangeCoeff = interParticleExchangeCoeff, iterPeriod = 10)


### set custom thickness to some range of sphere segments
updater.setThicknessToSpheres([0],0,3,0,7,0)# bIds, thetaMin, thetaMax, phiMin, phiMax, thickness
updater.setThicknessToSpheres([1],0,3,0,7,1)# bIds, thetaMin, thetaMax, phiMin, phiMax, thickness


######## engines
 
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb()]),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom6D()],
                [Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys(setCohesionOnNewContacts=True)],
                [Law2_ScGeom6D_CohFrictPhys_CohesionMoment()]
        ),
        NewtonIntegrator(gravity=(0, 0, 0), damping=0.0),
        updater,
]

O.dt = 1e-6
NstepsToRun = int(0.001/O.dt)+10


O.run(NstepsToRun, True)

if abs(O.bodies[0].state.coatingVolume[16] - 0.00024999999) > 1e-10: # the expected coating volume in given segment
	raise YadeCheckError("The volume of the particle's film is too far from expected.")
        
