if ('CGAL' in features):
	from yade import polyhedra_utils as p_util
	########## CONSTANTS (I add suffix '_poly' to variables to make them more unique (recommended for check scripts)
	edge_len_poly = 1  # dimension of the sample cuboid
	shift_poly = 5
	########## MATERIALS

	high_friction_id_poly = O.materials.append(PolyhedraMat(density=2500, frictionAngle=radians(89), label='high_friction_poly'))

	low_friction_id_poly = O.materials.append(PolyhedraMat(density=2500, frictionAngle=radians(1), label='low_friction_poly'))

	########## BODIES
	# Define base plates ( no. 1 with high friction and shifted no. 2 with low friction)
	bottom_plate_1_poly = p_util.polyhedra(
	        O.materials[high_friction_id_poly],
	        v=[
	                (-edge_len_poly, -edge_len_poly, -edge_len_poly / 2), (-edge_len_poly, edge_len_poly, -edge_len_poly / 2),
	                (-edge_len_poly, edge_len_poly, 0), (-edge_len_poly, -edge_len_poly, 0), (edge_len_poly, -edge_len_poly, -edge_len_poly / 2),
	                (edge_len_poly, edge_len_poly, -edge_len_poly / 2), (edge_len_poly, edge_len_poly, 0), (edge_len_poly, -edge_len_poly, 0)
	        ],
	        fixed=True
	)
	bottom_plate_2_poly = p_util.polyhedra(
	        O.materials[low_friction_id_poly],
	        v=[
	                (-edge_len_poly + shift_poly, -edge_len_poly, -edge_len_poly / 2), (-edge_len_poly + shift_poly, edge_len_poly, -edge_len_poly / 2),
	                (-edge_len_poly + shift_poly, edge_len_poly, 0), (-edge_len_poly + shift_poly, -edge_len_poly, 0),
	                (edge_len_poly + shift_poly, -edge_len_poly, -edge_len_poly / 2), (edge_len_poly + shift_poly, edge_len_poly, -edge_len_poly / 2),
	                (edge_len_poly + shift_poly, edge_len_poly, 0), (edge_len_poly + shift_poly, -edge_len_poly, 0)
	        ],
	        fixed=True
	)

	# Define samples  ( no. 1 with high friction and shifted no. 2 with low friction)
	sample_1_poly = p_util.polyhedra(
	        O.materials[high_friction_id_poly],
	        v=[
	                (-edge_len_poly / 2, -edge_len_poly / 2, 0), (-edge_len_poly / 2, edge_len_poly / 2, 0),
	                (-edge_len_poly / 2, edge_len_poly / 2, edge_len_poly), (-edge_len_poly / 2, -edge_len_poly / 2, edge_len_poly),
	                (edge_len_poly / 2, -edge_len_poly / 2, 0), (edge_len_poly / 2, edge_len_poly / 2, -0),
	                (edge_len_poly / 2, edge_len_poly / 2, edge_len_poly), (edge_len_poly / 2, -edge_len_poly / 2, edge_len_poly)
	        ],
	        fixed=False
	)
	sample_2_poly = p_util.polyhedra(
	        O.materials[low_friction_id_poly],
	        v=[
	                (-edge_len_poly / 2 + shift_poly, -edge_len_poly / 2, 0), (-edge_len_poly / 2 + shift_poly, edge_len_poly / 2, 0),
	                (-edge_len_poly / 2 + shift_poly, edge_len_poly / 2, edge_len_poly),
	                (-edge_len_poly / 2 + shift_poly, -edge_len_poly / 2, edge_len_poly), (edge_len_poly / 2 + shift_poly, -edge_len_poly / 2, 0),
	                (edge_len_poly / 2 + shift_poly, edge_len_poly / 2, -0), (edge_len_poly / 2 + shift_poly, edge_len_poly / 2, edge_len_poly),
	                (edge_len_poly / 2 + shift_poly, -edge_len_poly / 2, edge_len_poly)
	        ],
	        fixed=False
	)

	# Add bodies to simulation.
	bottom_plate_1_poly_id, bottom_plate_2_poly_id, sample_1_poly_id, sample_2_poly_id = O.bodies.append(
	        [bottom_plate_1_poly, bottom_plate_2_poly, sample_1_poly, sample_2_poly]
	)

	######### ENGINES
	O.engines = [
	        ForceResetter(),
	        InsertionSortCollider([Bo1_Polyhedra_Aabb()]),
	        InteractionLoop(
	                [Ig2_Polyhedra_Polyhedra_PolyhedraGeom()],
	                [Ip2_PolyhedraMat_PolyhedraMat_PolyhedraPhys()],
	                [Law2_PolyhedraGeom_PolyhedraPhys_Volumetric()],
	        ),
	        NewtonIntegrator(gravity=(0, 0, -10)),
	]

	######## set time step and run
	O.dt = 1e-4

	O.run(1000, wait=True)
	O.forces.setPermF(sample_1_poly_id, (1000, 0, 0))
	O.run(1000, wait=True)

	if O.interactions[bottom_plate_2_poly_id, sample_2_poly_id].phys.shearForce.norm() > 1e-5:
		raise YadeCheckError("Finished due to the 'leakage' in the shearForce between separate interactions!")
	else:
		print("OK - no 'leakage' in the shearForce between separate interactions.")
