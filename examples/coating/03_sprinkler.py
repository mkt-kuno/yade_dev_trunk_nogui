from yade import pack, export, timing

O.timingEnabled = True

# Visualisation disclaimer:
print("""\n\n Note that you will not be able to see the effect of sprinkler action directly in the Yade GUI. You will only see how the particles start agglomerate due to cohesion gained with thickness.
    To visualise the thickness of the exported files in ParaView:
    1) Apply programable filter copied from file programmable_filter_paraview.py
    2) Apply macro macro_ParaView_parametrized.py (remember about setting proper theta, phi, and range of visualised thickness.
    Alternatively you can use 'segSphExporter()' that is loaded from vtk_export_helpers.py, but the files generated are heavy and spheres are 'edgy' (converted to meshes).
    """)

### helper functions

def sumVolume():
    vol = 0
    for b in O.bodies:
        if isinstance(b.shape, Sphere):
            volumes = b.state.coatingVolume
            vol += sum(volumes)
    return vol
   
# loading helper functions from external file by just executing it
with open("vtk_export_helpers.py") as f:
    vtk_export_helpers = f.read()
exec(vtk_export_helpers) 

       
    
###
###### constatnts
nTheta = 4
nPhi = 8

exchangeMatWithFacets = False

### prepare exporter
home = os.path.expanduser("~")
output_dir = home +'/tmp/coating/output-sprinkler/'
if not(os.path.exists(output_dir )):
	os.makedirs(output_dir)
	
exporter = export.VTKExporter(output_dir )
export_strings = {'Quatx':'b.state.se3[1][0]','Quaty':'b.state.se3[1][1]','Quatz':'b.state.se3[1][2]','Quatw':'b.state.se3[1][3]'}

segSphExporter = SegmentedSphereExporter(filebase = output_dir, nTheta = nTheta, nPhi = nPhi , quality = 1)

for thetaNo in range(nTheta):
    for phiNo in range(nPhi):
        segmentPos = thetaNo * nPhi + phiNo
        export_strings.update({ "theta_{:d}_phi_{:d}".format(thetaNo,phiNo):"b.state.coatingThickness[{:d}]".format(segmentPos)})

def export_with_thickness():
    exporter.exportSpheres(what = export_strings)
    exporter.exportFacets(what = {"negSideCoating" : "b.state.coatingThickness[0]", "posSideCoating" : "b.state.coatingThickness[1]"})
    
####### BODIES
#### material
mat1 = CohFrictMatSeg(young = 1e7, secondaryYoung = 1e7,
                    normalCohesion = 0, shearCohesion = 0, secondaryNormalCohesion = 5e6, secondaryShearCohesion = 5e6, 
                    fragile = True,
                    minThickness = 0.1,
                    maxThickness = 3)#CohFrictMatSeg()

#### create rectangular box from facets
facets = O.bodies.append(geom.facetBox((.5, .5, .5), (.5, .5, .5), material = mat1))#, wallMask=31))

#### spheres
sp = pack.SpherePack()
sp.makeCloud((0, 0, 0), (1, 1, 1), rMean=.06, rRelFuzz=.01, num = 200)
spheres = sp.toSimulation(material = mat1)

###### ENGINES
Ip2 = Ip2_CohFrictMat_CohFrictMat_CohFrictPhys()
Ip2_seg = Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys(setCohesionOnNewContacts=True)

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()]),
        InteractionLoop(
                [Ig2_Facet_Sphere_ScGeom6D(), Ig2_Sphere_Sphere_ScGeom6D()],
                [Ip2_seg],#Ip2_seg
                [Law2_ScGeom6D_CohFrictPhys_CohesionMoment()]
        ),
        NewtonIntegrator(gravity=(0, 0, -10), damping=0.0),
        PyRunner(command='export_with_thickness()', virtPeriod=0.04, initRun = True),
        #PyRunner(command="segSphExporter()", virtPeriod=0.04, initRun = True), # If you have problems with visualisation, use this simple exporter (but output files will be heavy or low resolution).
        PyRunner(command='controller()', virtPeriod=0.1, initRun = False),
        RotationEngine(angularVelocity = 2, ids=facets, zeroPoint = (.5, .5, .5), rotationAxis = (0,1,0), rotateAroundZero = True, dead = True, label = 'rotEngine')
]


O.dt = .5 * PWaveTimeStep()


### INITIALIZE THICKNESS
# initialize updater
updater = SegmentedStateUpdater(thetaResolution = nTheta, phiResolution = nPhi, intraParticleExchangeCoeff = 500, interParticleExchangeCoeff = 5000, activateWettability = True,  iterPeriod = 10)# note that the material will net be exchanged periodically if iterPeriod or virtPeriod is not specified
# set thickness to spheres
updater.setThicknessToSpheres(spheres,0,3,0,7,0)# bIds, thetaMin, thetaMax, phiMin, phiMax, thickness
# initialize state facets state manually (initialization will be rather handled inside setter functions, but just testing)
for Id in facets: # Still - lack of initialization may cause some problems
    updater.initialize(Id,0)
    O.bodies[Id].state.allowMatExchange = exchangeMatWithFacets
    
    
# add updater to the engines 
O.engines += [updater]

### initialize sprinkler (SegmentedMatSprinkler())
sprinkler = SegmentedMatSprinkler(pos = (.5, .5, 2.0), alpha = 15, beta = 15, alphaResolution = 50, betaResolution = 50, feedRate = 100000, dead = True, virtPeriod = 5e-3)

O.engines += [sprinkler]

#### define controller

def controller():
    if sprinkler.dead and O.time > 2:
        sprinkler.dead = False
    if rotEngine.dead and O.time > 10:
        rotEngine.dead = False
        sprinkler.feedRate = 300000
    if O.time > 20:
        O.pause()
        timing.stats()
### run
O.run()
