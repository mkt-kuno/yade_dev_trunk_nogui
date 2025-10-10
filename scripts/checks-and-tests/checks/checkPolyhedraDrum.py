# 2025 © Karol Brzeziński<karol.brze @gmail.com>
#Rotating drum with polyhedra based on:
# 2020 © Vasileios Angelidakis<v.angelidakis2 @ncl.ac.uk>
# 2020 © Bruno Chareyre<bruno.chareyre @grenoble - inp.fr>
# 2020 © Robert Caulk<rob.caulk @gmail.com>

#Benchmark of basic performance of open - source DEM simulation systems
#Case 2 : Rotating Drum
"""
I added this check to detect 'Jumpy' particles.The particles in the rotating drum tend to achieve way to high velocity.                                        
Only happens for multithreading(-jN option where N > 1) with versions installed from binary files(everything is ok with versions compiled locally) .Hence, it is quite unlikely that this test will break CI pipeline. However, it will guard against such situation in the pipeline.
"""

if ('CGAL' in features):
	from yade import ymport, polyhedra_utils, pack
	import numpy as np

	#####################   0. HELPER FUNCTIONS  #####################


	def controller(startTime, checkVel=False):
		if rotation.dead and O.time > startTime:
			rotation.dead = False

		if checkVel:
			maxVel = 0
			for b in O.bodies:
				v = b.state.vel.norm()
				if v > maxVel:
					maxVel = v
#print("Maximum velocity of the particle is {:.1f} m/s".format(maxVel))
			if maxVel > 50:
				raise YadeCheckError("Finished due to the unstable ('jumpy') particles (velocity > 50 m/s)!")

	#####################   1. INPUT/OUTPUT  #####################

	angularVelocity = 2 * pi

	startTime = 0.0
	simulatedTime = 1.5

	#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- #
	#Materials

	Steel_poly = O.materials.append(PolyhedraMat(young=210e8, poisson=0.2, density=7200, label='Steel_poly'))

	M2_poly = O.materials.append(PolyhedraMat(young=0.5e8, poisson=0.2, density=2000, label='M2_poly'))

	#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- #

	#####################   2. YADE PART  #####################

	drum_width = 0.06
	drum_center = Vector3(0, 0, 0)
	drum_thickness = 0.005  # KB: assumed by myself
	drum_segments = 5
	drum_R_in = 0.1
	drum_R_out = drum_R_in + drum_thickness

	angles = np.linspace(0, 2 * pi, drum_segments)

	drum_polyhedra = []
	for angle in angles:
		vertices = []
		for y_side in [-1, 1]:
			for R in [drum_R_in, drum_R_out]:
				for angle_shift in [-pi / (drum_segments - 1), pi / (drum_segments - 1)]:
					y_part = Vector3(0, 0.5 * y_side * drum_width, 0)
					xz_part = R * Vector3(sin(angle + angle_shift), 0, cos(angle + angle_shift))
					vertices += [drum_center + y_part + xz_part]
		drum_polyhedra += [polyhedra_utils.polyhedra(O.materials[Steel_poly], v=vertices, color=(0, 1, 0), fixed=True)]
	drum_ids = O.bodies.append(drum_polyhedra)

	#polyhedra drum sides
	for side in [-1, 1]:
		vertices = [
		        Vector3(drum_R_out, side * drum_width / 2, drum_R_out),
		        Vector3(-drum_R_out, side * drum_width / 2, drum_R_out),
		        Vector3(drum_R_out, side * drum_width / 2, -drum_R_out),
		        Vector3(-drum_R_out, side * drum_width / 2, -drum_R_out),
		        Vector3(drum_R_out, side * (drum_width / 2 + drum_thickness), drum_R_out),
		        Vector3(-drum_R_out, side * (drum_width / 2 + drum_thickness), drum_R_out),
		        Vector3(drum_R_out, side * (drum_width / 2 + drum_thickness), -drum_R_out),
		        Vector3(-drum_R_out, side * (drum_width / 2 + drum_thickness), -drum_R_out)
		]
		drum_wall = polyhedra_utils.polyhedra(O.materials[Steel_poly], v=vertices, color=(0, 1, 0), fixed=True)
		drum_wall.shape.wire = True
		drum_ids += [O.bodies.append(drum_wall)]

	## Engines
	O.engines = [
	        ForceResetter(),
	        InsertionSortCollider([Bo1_Polyhedra_Aabb()]),
	        InteractionLoop(
	                [Ig2_Polyhedra_Polyhedra_PolyhedraGeom()],
	                [Ip2_PolyhedraMat_PolyhedraMat_PolyhedraPhys()],
	                [Law2_PolyhedraGeom_PolyhedraPhys_Volumetric()],
	        ),
	        NewtonIntegrator(damping=0.01, gravity=[0, 0, -9.81]),
	        RotationEngine(
	                dead=1,
	                rotateAroundZero=True,
	                zeroPoint=(0, 0, 0),
	                rotationAxis=(0, 1, 0),
	                angularVelocity=angularVelocity,
	                ids=drum_ids,
	                label='rotation'
	        ),
	        PyRunner(command="controller(startTime = startTime, checkVel = True)", virtPeriod=0.01, initRun=True)
	]

	#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- #

	#NSTEPS = 1000

	sp_temp = pack.SpherePack()
	mn_temp = -Vector3(drum_R_in, 0.5 * drum_width, drum_R_in)
	mx_temp = Vector3(drum_R_in, 0.5 * drum_width, drum_R_in)
	sp_temp.makeCloud(mn_temp, mx_temp, 0.002, 0, 200, False)

	cyl = pack.inCylinder((0, -0.5 * drum_width, 0), (0, 0.5 * drum_width, 0), drum_R_in * 0.5)
	sp_temp = pack.filterSpherePack(cyl, sp_temp, True)

	#replace spheres with polyhedra
	poly = []

	cube_vertices = [
	        Vector3(-1, -1, -1),
	        Vector3(1, -1, -1),
	        Vector3(-1, 1, -1),
	        Vector3(1, 1, -1),
	        Vector3(-1, -1, 1),
	        Vector3(1, -1, 1),
	        Vector3(-1, 1, 1),
	        Vector3(1, 1, 1)
	]

	##############
	for sph_temp in sp_temp:
		r_temp = sph_temp[1]
		pos_poly = sph_temp[0]
		mat_poly = O.materials[M2_poly]
		half_d = r_temp / 2 * (4 * pi / 3)**(1 / 3)
		v_poly = [vv_poly * half_d + pos_poly for vv_poly in cube_vertices]
		poly += [polyhedra_utils.polyhedra(material=mat_poly, v=v_poly)]

	O.bodies.append(poly)

	O.dt = 0.8 * PWaveTimeStep()

	O.run(int((simulatedTime + startTime) / O.dt), wait=True)

	print("OK - Drum with polyhedra finished without detecting 'jumpy' particles.")
