def phiFromStl(stlFile, nMin=80, extendCoeff=1):
	'''Return a (phi,grid) pair that can almost be passed to levelSetBody() (just think of turning phi numpy array into a list with .tolist()). Requires to "pip install trimesh rtree scipy". It is an alternative to using phiIniFromPts + FastMarchingMethod, which is much more computationally expensive (RAM requirements in particular may reach 10-100 of GB for a MB .stl and a usual 1e3 - 1e4 gridpoints grid) but certainly leads to less 0-values in LevelSet.distField
    :param string stlFile: full path to the .stl (file object also accepted ?). Has to conform inertial frame, it is possible to use mkStlInertial for this purpose
    :param nMin: minimum number of gridpoints discretizing the stl (just its mesh, whatever used *extendCoeff*) along an axis. The final grid will actually have a couple of more gridpoints than extendCoeff*nMin on that axis.
    :param extendCoeff: multiplicative coefficient to extend the grid beyond the mesh bounding box (if useful)
    '''
	try:
		numpy  # this one is for sure present alongside a YADE installation
	except NameError:  # we just need to be sure to have it imported
		import numpy
	try:
		import trimesh
	except ModuleNotFoundError:  # using except just for a print, not really avoiding the error
		print('****Python module trimesh is required, see the doc****')
		import trimesh  # we still repeat the same command for the error to appear to the user
	try:
		import rtree
	except ModuleNotFoundError:  # same philosophy as with trimesh
		print('****Python module rtree is required, see the doc****')
		import rtree
	mesh = trimesh.load_mesh(stlFile, process=False)  # process=False (as a trimesh.Trimesh attribute) for avoiding automatic merge of vertices
	aabb = mesh.bounding_box
	min_corner, max_corner = aabb.vertices.min(axis=0), aabb.vertices.max(axis=0)
	xmin, ymin, zmin = extendCoeff * min_corner
	xmax, ymax, zmax = extendCoeff * max_corner

	# Now copying-pasting the grid construction of phiIniFromPts (TODO: make a separate function and make the grid inertial ie with center of mass as a gridpoint)
	lxPt = xmax - xmin
	lyPt = ymax - ymin
	lzPt = zmax - zmin
	lMinPt = min(lxPt, lyPt, lzPt)

	# Deduced grid properties:
	spac = lMinPt / (nMin * extendCoeff)
	Nx, Ny, Nz = ceil(lxPt / spac), ceil(lyPt / spac), ceil(lzPt / spac)  # these are gridpoints strictly fitting the cloud
	nGP = Vector3i(Nx, Ny, Nz) + Vector3i(2, 2, 2)  # these will be the points of the actual RegularGrid, extended by one gp on each side:
	gridMin = Vector3(xmin, ymin, zmin) - spac * Vector3(1, 1, 1)
	grid = RegularGrid(min=gridMin, spacing=spac, nGP=nGP)

	# Preparing trimesh job:
	gpArray = numpy.array([[[grid.gridPoint(i, j, k) for k in range(nGP[2])] for j in range(nGP[1])] for i in range(nGP[0])])
	gpArrayReshaped = gpArray.reshape(
	        -1, 3
	)  # converting in correct format (with as many lines as necessary and 3 columns for x,y,z) as expected by trimesh below
	# Running trimesh to conclude:
	phiValues = -trimesh.proximity.signed_distance(
	        mesh, gpArrayReshaped
	)  # "-" for being positive outside, negative inside; with proximity.closest_point as the under-the-hood-instrumental function
	# phiValues is now a 1D-array with nGP values, we will convert its dimensions below
	return (phiValues.reshape(gpArray.shape[:-1]), grid)  # phiValues transformed to a 3D array (conforming gridpoints gpArray in terms of shape)
