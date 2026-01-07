from yade import plot, export

# Visualisation disclaimer:
print("""\n\n 
    To visualise the thickness of the exported files in ParaView:
    1) Apply programable filter copied from file programmable_filter_paraview.py
    2) Apply macro macro_ParaView_parametrized.py (remember about setting proper theta, phi, and range of visualised thickness.
    Alternatively you can use 'segSphExporter()' that is loaded from vtk_export_helpers.py, but the files generated are heavy and spheres are 'edgy' (converted to meshes).
    """)



### constants
nTheta = 8
nPhi = 16
interParticleExchangeCoeff = 1000000 # 
intraParticleExchangeCoeff = 100000 # 

### prepare exporter
home = os.path.expanduser("~")
output_dir = home +'/tmp/coating/output-example02/'
if not(os.path.exists(output_dir )):
	os.makedirs(output_dir)
	
exporter = export.VTKExporter(output_dir )
export_strings = {'Quatx':'b.state.se3[1][0]','Quaty':'b.state.se3[1][1]','Quatz':'b.state.se3[1][2]','Quatw':'b.state.se3[1][3]'}

for thetaNo in range(nTheta):
    for phiNo in range(nPhi):
        segmentPos = thetaNo * nPhi + phiNo
        export_strings.update({ "theta_{:d}_phi_{:d}".format(thetaNo,phiNo):"b.state.coatingThickness[{:d}]".format(segmentPos)})
        
        
### material
mat1 = CohFrictMatSeg()

### bodies
s1 = sphere((0,0,0), 1, material = mat1, fixed = True )
s2 = sphere((2.05,0,0), 1, material = mat1, fixed = True )
s1_id, s2_id = O.bodies.append([s1, s2])


### initialize updater
updater = SegmentedStateUpdater(thetaResolution = nTheta, phiResolution = nPhi, 
                                interParticleExchangeCoeff = interParticleExchangeCoeff, 
                                intraParticleExchangeCoeff = intraParticleExchangeCoeff,
                                activateWettability = True,
                                iterPeriod = 10)


### set custom thickness to some range of sphere segments
updater.setThicknessToSpheres([0],0,7,0,15,0)# bIds, thetaMin, thetaMax, phiMin, phiMax, thickness
updater.setThicknessToSpheres([1],7,7,0,15,1)# bIds, thetaMin, thetaMax, phiMin, phiMax, thickness


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
        PyRunner(command = 'exporter.exportSpheres(what = export_strings)', virtPeriod = 0.5, initRun = True),
        PyRunner(command = 'controller()', virtPeriod = 0.1),
]

O.dt = 1e-5
### functions
def controller():
    if O.time > 20 and s2.state.vel.norm() == 0:
        s2.state.vel = (-0.1,0,0)
    if O.time > 21 and s2.state.vel.norm() > 0:
        s2.state.vel = (0,0,0)
    if O.time >200:
        O.pause()

O.run()
        
