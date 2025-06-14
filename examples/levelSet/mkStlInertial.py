def mkStlInertial(stlFileName, verbose=False):
	'''Maps a stl mesh into one of its inertial frame. The new .stl will be saved, appending '_madeI.stl' to *stlFileName*. Requires stl Python module, see https://pypi.org/project/numpy-stl/.
    :stlFileName: path to the .stl to consider, wo the .stl extension (string)
    :verbose: boolean to control informative message outputs about the function process
    '''
	import numpy, time
	try:
		import stl  # = https://pypi.org/project/numpy-stl/. Alternate choice https://python-stl.readthedocs.io/en/latest/ is no longer maintained on https://github.com/apparentlymart/python-stl
	except ModuleNotFoundError:  # using except just for a print, not really avoiding the error
		print('****Python module stl is required, see the doc****')
		import stl  # we still repeat the same command for the error to appear to the user
	t1 = time.time()
	stlMesh = stl.mesh.Mesh.from_file(stlFileName + '.stl')
	t2 = time.time()
	if verbose:
		print('\nReading given file in', t2 - t1, 'seconds')
	volume, cog, inertia = stlMesh.get_mass_properties()  # cog being a numpy array
	t1 = time.time()
	if verbose:
		print('Computing mass properties in ~', t1 - t2, 'seconds')
		print('Initially, stl is centered at', cog, 'with inertia\n', inertia)

	# Getting eigenbasis information:
	[val, vect] = numpy.linalg.eig(inertia)  # vect is the eigenvectors basis in current basis
	if numpy.linalg.det(vect) < 0:  # eg (x,z,y)
		if verbose:
			print('Non right-handed basis, modifying vect, which is for now\n', vect)
		vect[:, 2] = -vect[:, 2]  # making it (x,z,-y) which is a more suitable basis
		if verbose:
			print('And now\n', vect, '\nwith determinant', numpy.linalg.det(vect))
	rot = vect.T  # taking the transpose that will be useful for changing basis
	# Applying the transformation, avoiding to use eg stlMesh.transform or stlMesh.rotate on purpose (with source code at https://github.com/WoLpH/numpy-stl/tree/develop/stl/base.py) because of https://github.com/WoLpH/numpy-stl/issues/166:
	for fIdx in range(len(stlMesh.vectors)):  # loop over facets
		for vIdx in range(3):  # loop over 3 vertices in each facet
			stlMesh.vectors[fIdx][vIdx] = rot @ (stlMesh.vectors[fIdx][vIdx] - cog)
	if verbose:
		volume, cog, inertia = stlMesh.get_mass_properties()
		print('\nAfter transformation, stl is centered at', cog, 'and has following inertia\n', inertia)
	stlMesh.save(stlFileName + '_madeI.stl')
