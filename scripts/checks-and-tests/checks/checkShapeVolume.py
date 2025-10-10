# Checking shape.getVolume()


# Handy function for comparing obtained quantities with expected ones (also used in a more rudimentary form in checkLS_DEM..)
def equalNbr(x, xRef, rtol=0.02):
	'''Whether x and xRef are equal numbers, up to a relative tolerance rtol.'''
	import numpy as np
	try:
		equal_bool = np.isclose(volume_obtained, volume_expected, rtol=rtol, atol=0)
	except TypeError:  # that would be the case in the ckHP* steps of the pipeline
		if abs(x - xRef) / abs(xRef) > rtol:
			equal_bool = False
		else:
			equal_bool = True
	return equal_bool


ctr = Vector3.Zero

### Box ###
box_extents = np.array([2, 1.3, 6.9])  # half-extents
bod = box(ctr, box_extents)
volume_expected = 8 * box_extents.prod()
volume_obtained = bod.shape.getVolume()
if not equalNbr(volume_obtained, volume_expected, rtol=1.e-8):
	raise YadeCheckError('Box.getVolume() gave', volume_obtained, 'vs', volume_expected, 'expected')

### Facet ###
bod = facet([ctr, Vector3.UnitX, Vector3.UnitY])
volume_obtained = bod.shape.getVolume()
if not (abs(volume_obtained) < 1.e-9):
	raise YadeCheckError('Facet.getVolume() gave', volume_obtained, 'vs 0 expected')

### LevelSet ###
if 'LS_DEM' in features:
	rad = 2.34
	bod = levelSetBody('sphere', ctr, rad)
	volume_expected = 53.692939593336376  # 53.67 for 4./3 * np.pi * rad**3
	volume_obtained = bod.shape.getVolume()
	if not equalNbr(volume_obtained, volume_expected, rtol=1.e-8):
		raise YadeCheckError('LevelSet.getVolume() gave', volume_obtained, 'vs', volume_expected, 'expected')

### Polyhedra ###
if 'CGAL' in features:
	import polyhedra_utils
	rad = 2.34
	bod = polyhedra_utils.polyhedralBall(rad, 20, FrictMat(), ctr)
	volume_expected = 39.27094028372987  # 53.67 for 4./3 * np.pi * rad**3
	volume_obtained = bod.shape.getVolume()
	if not equalNbr(volume_obtained, volume_expected, rtol=1.e-8):
		raise YadeCheckError('Polyhedra.getVolume() gave', volume_obtained, 'vs', volume_expected, 'expected')

### PotentialBlock ###
if 'POTENTIAL_BLOCKS' in features:
	import potential_utils
	bod = potential_utils.cuboid(FrictMat())  # 1x1x1 PB cuboid
	volume_expected = 0.9999999999999998
	volume_obtained = bod.shape.getVolume()
	if not equalNbr(volume_obtained, volume_expected, rtol=1.e-8):
		raise YadeCheckError('PotentialBlock.getVolume() gave', volume_obtained, 'vs', volume_expected, 'expected')

### Sphere ###
rad = 2.34
bod = sphere(ctr, rad)
volume_expected = 4. / 3 * np.pi * rad**3
volume_obtained = bod.shape.getVolume()
if not equalNbr(volume_obtained, volume_expected, rtol=1.e-8):
	raise YadeCheckError('Sphere.getVolume() gave', volume_obtained, 'vs', volume_expected, 'expected')

### Wall ###
bod = wall(0, 1)
volume_obtained = bod.shape.getVolume()
if not (abs(volume_obtained) < 1.e-9):
	raise YadeCheckError('Wall.getVolume() gave', volume_obtained, 'vs 0 expected')
