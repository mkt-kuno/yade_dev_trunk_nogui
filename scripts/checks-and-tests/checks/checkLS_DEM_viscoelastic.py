# -*- encoding=utf-8 -*-
# Vasileios Angelidakis <v.angelidakis@qub.ac.uk>

if ('LS_DEM' in features):

	#---------------------------------------------------------------------
	# Function to compare if two numerical values are nearly equal, within a given relative tolerance
	def equalNbr(x, xRef, tol=0.001):
		'''Whether x and xRef are equal numbers, up to a relative tolerance tol'''
		if abs(x - xRef) / abs(xRef) > tol:
			return False
		else:
			return True

	#---------------------------------------------------------------------
	# Material
	O.materials.append(ViscElMat(young=-1, kn=4000, ks=4000, en=0.6, et=0.6, poisson=0.4, density=1000, frictionAngle=atan(0.5), label='particles'))

	#---------------------------------------------------------------------
	# Free body and fixed body
	O.bodies.append(levelSetBody('sphere', (0, 0, 0.03), 0.01, spacing=0.0005, material='particles', nodesPath=1, nSurfNodes=83))
	O.bodies.append(levelSetBody('sphere', (0, 0, 0), 0.01, spacing=0.0005, material='particles', nodesPath=1, nSurfNodes=83))
	O.bodies[-1].state.blockedDOFs = 'xyzXYZ'

	#---------------------------------------------------------------------
	# Engines
	O.engines = [
	        ForceResetter(),
	        InsertionSortCollider([Bo1_LevelSet_Aabb()], verletDist=0.0, label="collider"),
	        InteractionLoop(
	                [Ig2_LevelSet_LevelSet_MultiScGeom(label='ig2')], [Ip2_ViscElMat_ViscElMat_MultiViscElPhys(label='ip2')],
	                [Law2_MultiScGeom_MultiViscElPhys_Basic(label='law')],
	                label='ILoop'
	        ),
	        NewtonIntegrator(damping=0.0, gravity=(0, 0, 0), label='newton'),
	]

	O.dt = 1e-6

	#---------------------------------------------------------------------
	# Set initial velocity before collision
	b = O.bodies[0]
	v0 = -0.4
	b.state.vel = [0, 0, v0]

	#---------------------------------------------------------------------
	# Run enough timesteps for a contact to happen, and compare the requested coefficient of restitution with the observed one.
	O.run(50000, True)
	v1 = O.bodies[0].state.vel[2]

	if not equalNbr(abs(v1 / v0), O.materials[0].en, 0.001):
		raise YadeCheckError(
		        "Failed because the observed coefficient of restitution:", abs(v1 / v0), "is not the same as the expected value of", O.materials[0].en
		)

else:
	print("Skip checkLS_DEM_viscoelastic, LS-DEM feature not available")
