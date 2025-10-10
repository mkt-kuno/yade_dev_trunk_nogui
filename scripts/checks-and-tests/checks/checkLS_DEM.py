if ('LS_DEM' in features):
	# Handy functions for comparing obtained quantities with expected ones
	######################################################################
	def equalVectors(n, nRef, tol=0.02):
		'''Whether n and nRef are equal or opposite vectors, up to a relative tolerance tol (that applies to the norm)'''
		if (n - nRef).norm() > tol * nRef.norm() and (n + nRef).norm() > tol * nRef.norm():
			return False
		else:
			return True

	def equalNbr(x, xRef, tol=0.02):
		'''Whether x and xRef are equal numbers, up to a relative tolerance tol'''
		if abs(x - xRef) / abs(xRef) > tol:
			return False
		else:
			return True

	# Particle-scale comparisons for a sphere
	#########################################

	lsSph = levelSetBody('sphere', radius=1, spacing=0.05, nodesPath=1, smearCoeff=1)

	voxExpected = 33371  # 33371 grid voxels being inside (as detected by below numpy.sum operation) at the introduction of .binarize()
	voxObtained = numpy.sum(numpy.array(lsSph.shape.binarize()))
	if voxObtained != voxExpected:
		raise YadeCheckError("Incorrect LevelSet.binarize() with", voxObtained, "voxels being detected as inside, vs", voxExpected, "expected")

	def distSphereTh(pt, radius=1):
		return Vector3(pt).norm() - radius

	pts = [Vector3(0.11, 0.02, 0.08), 2 * Vector3.Ones, 100 * Vector3.Ones]
	distVals = [distSphereTh(pt) for pt in pts]
	normals = [pt / pt.norm() for pt in pts]  # Normal of a sphere is c*(x,y,z) where c is the normalisation constant.
	if not equalNbr(lsSph.shape.distance(pts[0]), distVals[0], 5.e-3):
		raise YadeCheckError(
		        "Incorrect distance value to a unit sphere for an inside point", pts[0], ":", lsSph.shape.distance(pts[0]), "vs", distVals[0],
		        "expected."
		)
	if not equalNbr(lsSph.shape.distance(pts[0], True), distVals[0], 5.e-3):
		raise YadeCheckError(
		        "When using unbound=True, incorrect distance value to a unit sphere for an inside point", pts[0], ":", lsSph.shape.distance(pts[0]),
		        "vs", distVals[0], "expected."
		)
	if not equalVectors(lsSph.shape.normal(pts[0]), normals[0], 0.06):
		raise YadeCheckError(
		        "Incorrect normal for a unit sphere for an inside point", pts[0], ":", lsSph.shape.normal(pts[0]), "vs", normals[0], "expected."
		)
	for idx in [1, 2]:  # out-of-the grid points
		dist = lsSph.shape.distance(pts[idx], True)
		distTh = distVals[idx]
		if not equalNbr(dist, distTh, 5.e-3):
			raise YadeCheckError("Incorrect distance value to a unit sphere for an outside point", pts[idx], ":", dist, "vs", distTh, "expected.")
		normalCalc = lsSph.shape.normal(pts[idx], True)
		normalTh = normals[idx]
		if not equalVectors(normalCalc, normalTh, 0.06):
			raise YadeCheckError("Incorrect normal for a unit sphere for an outside point", pts[idx], ":", normalCalc, "vs", normalTh, "expected.")
	volTh = 4. / 3. * pi  # expected volume
	if not equalNbr(lsSph.shape.getVolume(), volTh, 4.e-3):
		raise YadeCheckError("Failed because of an incorrect sphere volume in LS-DEM:", lsSph.shape.getVolume(), "vs", volTh, "expected")
	lNorm = [nod.norm() for nod in lsSph.shape.surfNodes]
	if not equalNbr(max(lNorm), 1, 1.e-3) or not equalNbr(min(lNorm), 1, 1.e-3):
		raise YadeCheckError("Failed because of incorrect boundary nodes on a sphere in LS-DEM")
	surf, surfTh = lsSph.shape.getSurface(), 4 * pi
	if not equalNbr(surf, surfTh, 1.2e-3):
		raise YadeCheckError("Incorrect getSurface() for a sphere in LS-DEM: got ", surf, "vs", surfTh)

	# Particle-scale comparisons for a superellipsoid
	#################################################

	### With respect to exact (volume) result:
	rx, ry, rz, epsE, epsN = 0.5, 1.2, 1., 0.1, 0.5
	lsSe = levelSetBody('superellipsoid', extents=(rx, ry, rz), epsilons=(epsE, epsN), spacing=0.05)

	# see A.H. Barr, in Graphics Gems III, D. Kirk (1995) for the following expressions of a superellipsoid volume. With log therein = natural logarithm
	def funG(x):  # see also http://people.math.sfu.ca/~cbm/aands/abramowitz_and_stegun.pdf p. 75
		if x <= 0:
			raise YadeCheckError("Gamma function does not apply to", x, "< 0")
		gam0, gam1, gam2, gam3 = 1. / 12, 1. / 30, 53. / 210, 195. / 371
		gam4, gam5 = 22999. / 22737, 29944523. / 19733142
		gam6 = 109535241009. / 48264275462
		contFrac = gam0 / (x + gam1 / (x + gam2 / (x + gam3 / (x + gam4 / (x + gam5 / (x + gam6 / x))))))
		return 0.5 * log(2 * pi) - x + (x - 0.5) * log(x) + contFrac

	def funGamma(x):  # exp(funG) is already the gamma function. Using nevertheless the factorial property for a more precise evaluation as below
		return exp(funG(x + 5)) / (x * (x + 1) * (x + 2) * (x + 3) * (x + 4))

	def beta(x, y):
		return funGamma(x) * funGamma(y) / funGamma(x + y)

	volExp = 2. / 3. * rx * ry**rz * epsE * epsN * beta(epsE / 2, epsE / 2) * beta(epsN, epsN / 2)
	if not equalNbr(lsSe.shape.getVolume(), volExp, 0.05):
		raise YadeCheckError("Failed because of an incorrect superellipsoid volume in LS-DEM:", lsSe.shape.getVolume(), "vs", volExp, "expected")

	### With respect to previous YADE-obtained results (this is ~ the script of https://gitlab.com/yade-dev/trunk/-/issues/375):
	rx, ry, rz, epsE, epsN = [0.4, 1., 0.8, 0.4, 1.6]  # Shape E from Duriez2021b = Duriez & Galusinski (2021) Computers & Geosciences 157
	volTh, inertiaTh = 1.0864026569757073, numpy.array(
	        [0.318409979166153656, 0.1283513094313180336, 0.2624619724909635354]
	)  # theoretical volume and inertia (xx,yy,zz) coefficients/eigenvalues as per Barr1995 and reported by Duriez 2021b in Table 2. Mind the typo therein which inverted Ixx and Iyy for that particular shape..
	#### without smearing, with reference error data from Duriez2021b (Figs. 7 and 15):
	resVals = [4, 5, 8]  # considered values for grid resolution
	volRefError = numpy.array([0.97937905, 1.02173535, 1.0014703])  # for grid resolution 4, 5 and 8
	inertiaRefError = numpy.array(
	        [
	                [1.05002943, 0.98222898, 1.00315571]  # for res = 4
	                ,
	                [1.02746689, 1.03917029, 1.03074973]  # for res = 5
	                ,
	                [1.00153711, 1.0013097, 0.99237233]  # for res = 8
	        ]
	)
	# Preparing the array for presently obtained errors
	volObtainedError = -numpy.ones(len(resVals))  # 1 item of volume error per grid resolution
	inertiaObtainedError = -numpy.ones((len(resVals), 3))  # 3 items of inertia errors per grid resolution
	# the LS description:
	for idx, res in enumerate(resVals):
		b = levelSetBody("superellipsoid", extents=(rx, ry, rz), epsilons=(epsE, epsN), spacing=2 * min(rx, ry, rz) / res, nSurfNodes=0, smearCoeff=-1)
		volLS, inertiaLS = b.shape.getVolume(), numpy.array(b.shape.inertia())
		volObtainedError[idx] = volLS / volTh
		inertiaObtainedError[idx, :] = inertiaLS / inertiaTh
	# the comparisons, with a higher tolerance for the finer grid to accept numeric variation in a FAST_NATIVE build:
	if not numpy.all(numpy.isclose(volObtainedError[:2], volRefError[:2], rtol=1.e-5)
	                ) or not numpy.isclose(volObtainedError[2], volRefError[2], rtol=1.e-3):
		raise YadeCheckError(
		        'Problem on LS volume (inertia not yet tested) of Duriez2021b shape E, got following ratios for the different grid resolutions\n',
		        volObtainedError, 'vs\n', volRefError, 'expected, wo smearing'
		)
	if not numpy.all(numpy.isclose(inertiaObtainedError[:2], inertiaRefError[:2], rtol=1.e-5)
	                ) or not numpy.all(numpy.isclose(inertiaObtainedError[:2], inertiaRefError[:2], rtol=5.e-3)):
		raise YadeCheckError(
		        'Problem on LS inertia (volume is OK) of Duriez2021b shape E, got following ratios for the different grid resolutions\n',
		        inertiaObtainedError, 'vs\n', inertiaRefError, 'expected, wo smearing'
		)
	#### with smearing:
	b = levelSetBody("superellipsoid", extents=(rx, ry, rz), epsilons=(epsE, epsN), spacing=2 * min(rx, ry, rz) / 5, nSurfNodes=0, smearCoeff=1)
	# reference values from Yade 2025-05-13.git-5018262:
	volExpected, inertiaExpected = 1.1241686219537839, numpy.array(Vector3(0.3439694527863318974, 0.1415271320033010538, 0.2850711020402842966))
	# LS values and comparisons:
	volLS, inertiaLS = b.shape.getVolume(), numpy.array(b.shape.inertia())
	if not numpy.isclose(volLS, volExpected):
		raise YadeCheckError('Problem on LS volume (inertia not yet tested) of Duriez2021b shape E with smearing, got', volLS, 'vs', volExpected)
	if not numpy.all(numpy.isclose(inertiaLS, inertiaExpected)):
		raise YadeCheckError('Problem on LS inertia of Duriez2021b shape E with smearing, got', inertiaLS, 'vs', inertiaExpected)

	print('LS-DEM as correct as expected at particle scale')

	# Now looking at the relative movements of 2 spheres and 2 LevelSet-shaped twins
	################################################################################

	rad = 1  # the smallest sphere
	rRatio = 1.8  # rBig / rSmall
	centrSmall, centrBig = (0, 0, 0), (0, 0, rad * (1 + rRatio))
	prec = 80  # grid fineness
	nSurfNodes = 2502
	# the 2 true spheres, along z-axis
	O.bodies.append(sphere(centrSmall, rad, dynamic=False))
	O.bodies.append(sphere(centrBig, rRatio * rad, dynamic=False))
	movSph = O.bodies[1]

	# their LevelSet counterparts:
	O.bodies.append(levelSetBody('sphere', centrSmall, rad, spacing=2 * rad / prec, nSurfNodes=nSurfNodes, nodesPath=1, dynamic=False))
	O.bodies.append(levelSetBody('sphere', centrBig, rRatio * rad, spacing=2 * rRatio * rad / prec, nSurfNodes=nSurfNodes, nodesPath=1, dynamic=False))
	movLS = O.bodies[3]

	O.engines = [
	        ForceResetter(),
	        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_LevelSet_Aabb()]),
	        InteractionLoop(
	                [Ig2_Sphere_Sphere_ScGeom(avoidGranularRatcheting=False),
	                 Ig2_LevelSet_LevelSet_ScGeom()],
	                [Ip2_FrictMat_FrictMat_FrictPhys(kn=MatchMaker(algo='val', val=1.e7), ks=MatchMaker(algo='val', val=1.e7))],
	                [Law2_ScGeom_FrictPhys_CundallStrack(sphericalBodies=False)]
	        ),
	        NewtonIntegrator()
	]
	O.dt = 5.e-4

	# 1. Normal relative displacement
	movSph.state.vel = movLS.state.vel = (0, 0, -1.)
	O.run(100, True)

	if not O.interactions[2, 3]:
		raise YadeCheckError("No LS interaction after first stage of relative normal displacement")
	lsCont = O.interactions[2, 3]
	sphCont = O.interactions[0, 1]

	if not equalNbr(lsCont.geom.penetrationDepth, sphCont.geom.penetrationDepth, 1.e-12):  # 2.2e-13 is a feasible goal on that ideal case
		raise YadeCheckError(
		        "Normal overlap is too wrong in LS-DEM after first stage:", lsCont.geom.penetrationDepth, "vs", sphCont.geom.penetrationDepth, "in DEM"
		)

	if not equalVectors(lsCont.geom.normal, sphCont.geom.normal):
		raise YadeCheckError("Normal vector excessive mismatch after first stage:", sphCont.geom.normal, "vs", lsCont.geom.normal)
	movSph.state.vel = movLS.state.vel = Vector3.Zero

	#2. Circular relative displacement (~ pure shear)
	arc = rad * (1 + rRatio)
	dAlpha = 1.2 * pi / 4.
	nIt = 1500
	normV = arc * dAlpha / (nIt * O.dt)
	lsShearDisp = sphShearDisp = Vector3.Zero
	for i in range(nIt):
		movSph.state.vel = movLS.state.vel = normV * Vector3(cos(dAlpha / nIt * i), 0, -sin(dAlpha / nIt * i))
		O.step()
		lsShearDisp += lsCont.geom.shearInc
		sphShearDisp += sphCont.geom.shearInc

	if not equalNbr(
	        lsCont.geom.penetrationDepth, sphCont.geom.penetrationDepth, 0.03
	):  # 0.0284 error is expected here, would be eg 0.007 with 6402 nodes and grid precision 80
		raise YadeCheckError("Normal overlaps are too different after 2nd stage:", lsCont.geom.penetrationDepth, "vs", sphCont.geom.penetrationDepth)
	if not equalVectors(lsShearDisp, sphShearDisp):
		raise YadeCheckError("Excessive error on shear displacement after second stage:", sphShearDisp, "vs", lsShearDisp)
	if not equalVectors(
	        lsCont.geom.normal, sphCont.geom.normal, 0.03
	):  # allowing here 3 % of error. 6402 nodes and grid precision 80 would allow to go under 2 %
		raise YadeCheckError("Excessive normal vector mismatch after second stage:", sphCont.geom.normal, "vs", lsCont.geom.normal)
	print('LS-DEM ScGeom contact description as correct as expected')

	# A similar test about 2 LevelSet-spheres in relative movement, in a LSnodeGeom fashion and wo the truth twins
	##############################################################################################################

	# (jduriez note: this is zi2ter)

	O.reset()
	O.bodies.append(levelSetBody('sphere', centrSmall, rad, spacing=2 * rad / prec, nSurfNodes=nSurfNodes, nodesPath=2, n_neighborsNodes=6, dynamic=False))
	O.bodies.append(
	        levelSetBody(
	                'sphere',
	                centrBig,
	                rRatio * rad,
	                spacing=2 * rRatio * rad / prec,
	                nSurfNodes=nSurfNodes,
	                n_neighborsNodes=6,
	                nodesPath=2,
	                dynamic=False
	        )
	)
	movLS = O.bodies[1]

	O.engines = [
	        ForceResetter(),
	        InsertionSortCollider([Bo1_LevelSet_Aabb()], verletDist=0),
	        InteractionLoop(
	                [Ig2_LevelSet_LevelSet_LSnodeGeom(label='ig2')],
	                [Ip2_FrictMat_FrictMat_FrictPhys(kn=MatchMaker(algo='val', val=1.e7), ks=MatchMaker(algo='val', val=1.e7))],
	                [Law2_ScGeom_FrictPhys_CundallStrack(sphericalBodies=False)]
	        ),
	        NewtonIntegrator()
	]
	O.dt = 5.e-4

	# 1. Normal relative displacement
	movLS.state.vel = (0, 0, -1.)
	O.run(100, True)

	if not O.interactions[0, 1]:
		raise YadeCheckError("No LS(nodeGeom) interaction after first stage of relative normal displacement")
	lsCont = O.interactions[0, 1]
	unRef, nodeIdxRef = 0.04855479185286475, nSurfNodes - 1  # NB: we changed path since the previous part of the test so expected values have changed as well

	nMismatchExpected = 0.05587214056599649
	nMismatch = ig2.normalMismatch(lsCont)
	if not equalNbr(nMismatch, nMismatchExpected):
		raise YadeCheckError("ig2.normalMismatch gives", nMismatch, "vs", nMismatchExpected, "expected")

	if not equalNbr(lsCont.geom.penetrationDepth, unRef, 1.e-12):
		raise YadeCheckError("Normal overlap with LSnodeGeom after first stage is:", lsCont.geom.penetrationDepth, "vs", unRef, "expected")
	if not (lsCont.geom.surfNodeIdx == nodeIdxRef):
		raise YadeCheckError(
		        "Contacting node after first stage was detected as", lsCont.geom.surfNodeIdx, "vs", nodeIdxRef, "expected (about the pole)"
		)

	#2. Circular relative displacement (~ pure shear)
	arc = rad * (1 + rRatio)
	dAlpha = 1.2 * pi / 4.
	nIt = 1500
	normV = arc * dAlpha / (nIt * O.dt)
	lsShearDisp = sphShearDisp = Vector3.Zero
	for i in range(nIt):
		movLS.state.vel = normV * Vector3(cos(dAlpha / nIt * i), 0, -sin(dAlpha / nIt * i))
		O.step()
		lsShearDisp += lsCont.geom.shearInc
	shearDispRef, nodeIdxRef = Vector3(2.251877401227810882, -8.274556988554289828e-05, -1.172519769676061863), 1974
	if not equalVectors(lsShearDisp, shearDispRef):
		raise YadeCheckError("Obtained shear displacement with LSnodeGeom is:", lsShearDisp, "vs", shearDispRef)
	if not lsCont.geom.surfNodeIdx == nodeIdxRef:
		raise YadeCheckError("Contacting node after second stage is", lsCont.geom.surfNodeIdx, "vs", nodeIdxRef, "expected")

	print('LS-DEM LSnodeGeom contact description as correct as expected')

	# Wall - LevelSet test
	######################

	# jduriez note: this is lsc2
	O.reset()
	O.bodies.append(wall(0, 2))
	O.bodies.append(
	        levelSetBody(
	                'superellipsoid', (0, 0, 1.9),
	                extents=(2, 1, 1),
	                epsilons=(1, 1),
	                spacing=1. / 10,
	                nSurfNodes=3602,
	                nodesPath=1,
	                orientation=Quaternion((0, 1, 0), pi / 2.),
	                dynamic=False
	        )
	)
	O.engines = [
	        ForceResetter(),
	        InsertionSortCollider(
	                [Bo1_LevelSet_Aabb(), Bo1_Wall_Aabb()], verletDist=0
	        )  # no need to wait for IScollider to decide himself to use 0. Note that I can not use anything else for now: the bounds have to match the grids, and those do not extend much more than the bodies...
	        ,
	        InteractionLoop([Ig2_Wall_LevelSet_ScGeom()], [Ip2_FrictMat_FrictMat_FrictPhys()], [Law2_ScGeom_FrictPhys_CundallStrack(label='csLaw2')])
	]
	O.dt = 4.e-5
	O.step()
	if not O.interactions.has(0, 1, True):
		raise YadeCheckError('Expected Wall-LevelSet interaction not found')
	else:
		c_geom = O.interactions[0, 1].geom  # just geom is already a module name
		if not equalNbr(c_geom.penetrationDepth, 0.1, 0.051):
			raise YadeCheckError('Overlap not correctly measured here: got', c_geom.penetrationDepth, 'vs 0.1 expected with a 5.1% tolerance')
		O.bodies[1].state.pos = Vector3(0, 0, 0.5)  # not recommended elsewhere
		O.step()
		if not equalNbr(c_geom.penetrationDepth, 1.5, 0.0041):
			raise YadeCheckError('Overlap not correctly measured there: got', c_geom.penetrationDepth, 'vs 1.5 expected with a 0.41% tolerance')
		O.bodies[1].state.pos = Vector3(1, -3, -1.9)
		O.step()
		if not equalNbr(c_geom.penetrationDepth, 0.1, 0.052):
			raise YadeCheckError('Overlap not correctly measured there: got', c_geom.penetrationDepth, 'vs 0.1 expected with a 5.2% tolerance')
	print('Ig2_Wall_LevelSet_ScGeom as correct as expected')

	# Fast Marching Method applied to the distance to the unit sphere
	#################################################################
	grid = RegularGrid(-1.1, 1.1, 23)  # a cubic grid from -1.1 to 1.1 with 23 gp ie a 0.1 step
	errorExpected = 0.009393853398395624  # e.g. on Ubuntu 20.04.3 and jduriez axp17* while would be 0.009419100794945902 on Debian Bullseye, see https://gitlab.com/yade-dev/trunk/-/jobs/1832563583
	for heapChoice in [True, False]:
		fmm = FastMarchingMethod(heapSort=heapChoice, phiIni=distIniSE(radii=[1, 1, 1], epsilons=[1, 1], grid=grid), grid=grid)
		phiField = fmm.phi()
		error = 0
		for i in range(23):
			for j in range(23):
				for k in range(23):
					phi_ijk = phiField[i][j][k]
					error += abs(grid.gridPoint(i, j, k).norm() - 1 - phi_ijk)
		error /= grid.nGP.prod()  # average (and dimensionless with respect to unit radius)
		if not equalNbr(error, errorExpected, 5.e-3):
			raise YadeCheckError("Fast Marching Method with heapSort =", heapChoice, "gives an error of", error, "vs", errorExpected, "expected")
	print('Fast Marching Method as correct as expected')

	# Checking LS-twins of (multi-sphere) clumps
	############################################

	O.reset()
	O.bodies.appendClumped([sphere((0, 0, 0), 1)] + [sphere((2, 0, -4 + k * 2), 1) for k in range(5)])
	O.bodies.append(levelSetBody(clump=O.bodies[6].shape, spacing=0.1, nSurfNodes=0))
	vLSclump = O.bodies[7].shape.getVolume()
	vclump = 6 * 4. / 3 * pi  # = the one of the Clump body, indeed
	if not equalNbr(vLSclump, vclump, 5.e-4):
		raise YadeCheckError("LStwin of a 6 unit-spheres clump has a volume of", vLSclump, "vs", vclump, "expected")
	print('LStwin of clumps as correct as expected')
else:
	print("Skip checkLSdem, LS-DEM feature not available")
