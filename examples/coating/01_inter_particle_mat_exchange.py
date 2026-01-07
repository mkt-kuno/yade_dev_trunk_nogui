from yade import plot

# Visualisation disclaimer:
print("""\n\n 
    To visualise the thickness of the exported files in ParaView:
    1) Apply programable filter copied from file programmable_filter_paraview.py
    2) Apply macro macro_ParaView_parametrized.py (remember about setting proper theta, phi, and range of visualised thickness.
    Alternatively you can use 'segSphExporter()' that is loaded from vtk_export_helpers.py, but the files generated are heavy and spheres are 'edgy' (converted to meshes).
    """)


### constants
nTheta = 4
nPhi = 8
interParticleExchangeCoeff = 25
plotVirtPeriod = 0.01 # 

### material
mat1 = CohFrictMatSeg()

### bodies
s1 = sphere((0,0,0), 1, material = mat1, fixed = True )
s2 = sphere((2.05,0,0), 1, material = mat1, fixed = True )
s1_id, s2_id = O.bodies.append([s1, s2])

s2.state.vel = (-0.1,0,0)

### initialize updater
updater = SegmentedStateUpdater(thetaResolution = nTheta, phiResolution = nPhi, interParticleExchangeCoeff = interParticleExchangeCoeff, iterPeriod = 10)


### set custom thickness to some range of sphere segments
updater.setThicknessToSpheres([0],0,3,0,7,0)# bIds, thetaMin, thetaMax, phiMin, phiMax, thickness
updater.setThicknessToSpheres([1],0,3,0,7,1)# bIds, thetaMin, thetaMax, phiMin, phiMax, thickness


######## engines
Ig2 = Ig2_Sphere_Sphere_ScGeom6D()
Ip2_seg = Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys(setCohesionOnNewContacts=True)

Law2 = Law2_ScGeom6D_CohFrictPhys_CohesionMoment()
    
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb(),Bo1_Box_Aabb()]),
        InteractionLoop(
                [Ig2],
                [Ip2_seg],
                [Law2]
        ),
        NewtonIntegrator(gravity=(0, 0, 0), damping=0.0),
        updater,
        PyRunner(command = 'controller()', virtPeriod = 1),
        PyRunner(command = 'addPlotData()', virtPeriod = plotVirtPeriod)
]

O.dt = 1e-6
### functions
def controller():
    # revert velocity of sphere2
    s2.state.vel *= -1
    if O.time > 8:
        O.pause()
        
def addPlotData():
    vol0 = O.bodies[0].state.coatingVolume[16] #I just checked number of contacting segments manually
    vol1 = O.bodies[1].state.coatingVolume[20]
    if len(plot.data['t']) > 0:
        vol0_change_rate = (vol0-plot.data["vol0"][-1]) / plotVirtPeriod
        vol1_change_rate = (vol1-plot.data["vol1"][-1]) / plotVirtPeriod
        vol0_tot_change = vol0-plot.data["vol0"][0] 
        vol1_tot_change = vol1-plot.data["vol1"][0] 
    else:
        vol0_change_rate = 0
        vol1_change_rate = 0
        vol0_tot_change = 0
        vol1_tot_change = 0
        
    t = O.time
    
    plot.addData(t = t, t_ = t, vol0 = vol0, vol1 = vol1, vol0_change_rate = vol0_change_rate, vol1_change_rate = vol1_change_rate, vol0_tot_change = vol0_tot_change, vol1_tot_change = vol1_tot_change)
    
plot.plots = {'t':('vol0_change_rate', 'vol1_change_rate'), 't_':('vol1_tot_change')}
plot.plot(subPlots = True)

O.run()
        
