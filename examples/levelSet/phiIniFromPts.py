def phiIniFromPts(ptsX, ptsY, ptsZ, nMin=80, sclFactor=1):
	# Initial version was developed by cedric.galusinski@univ-tln.fr
	'''Return a (phiIni,grid) pair for an appropriate initialization of a Fast Marching Method to run in order to eventually construct a distance to a surface defined as a cloud of points. It can act as an alternative to using phiFromStl, the present workflow being much less computationally expensive in terms of RAM and time and does not involve additional external requirements, at the cost of having maybe too many 0-values in LevelSet.distField
    :param ptsX: 1D iterable (typically list, tuple, or numpy array) with the x-coordinates of all points
    :param ptsY: 1D iterable with the y-coordinates of all points
    :param ptsZ: 1D iterable with the z-coordinates of all points
    :param nMin: minimum number of gridpoints discretizing the cloud along an axis. The final grid will have a couple of more gridpoints. Algorithm failures may happen if nMin is too big with respect to the resolution of given cloud.
    :param sclFactor: a possible multiplying factor to apply on given coordinates, e.g. to adapt the unit system
    :return: a (phiIni,grid) tuple with phiIni as a 3D numpy array full of -inf, inf or 0, to turn into a list with .tolist() and grid as RegularGrid instance
    >>> r=2 # theoretical sphere volume is 33.510321638291124
    >>> ptsX=[1+r*sin(radians(theta))*cos(radians(phi)) for theta in range(0,360,1) for phi in range(0,180,1)]
    >>> ptsY=[2+r*sin(radians(theta))*sin(radians(phi)) for theta in range(0,360,1) for phi in range(0,180,1)]
    >>> ptsZ=[1+r*cos(radians(theta)) for theta in range(0,360,1) for phi in range(0,180,1)]
    >>> (phiIni,grid) = phiIniFromPts(ptsX,ptsY,ptsZ,nMin=60)
    >>> numpy.sum(phiIni<0)*grid.spacing**3
    31.482074074074074
    >>> numpy.sum(phiIni<=0)*grid.spacing**3
    35.517037037037035
    '''
	try:  # following tests should directly work in YADE
		numpy
		math
	except NameError:
		import numpy, math
	ptsX = numpy.fromiter(ptsX, 'double') * sclFactor
	ptsY = numpy.fromiter(ptsY, 'double') * sclFactor
	ptsZ = numpy.fromiter(ptsZ, 'double') * sclFactor

	# Point cloud properties in terms of Axis Aligned Bounding Box:
	xmin = min(ptsX)
	ymin = min(ptsY)
	zmin = min(ptsZ)
	xmax = max(ptsX)
	ymax = max(ptsY)
	zmax = max(ptsZ)
	lxPt = xmax - xmin
	lyPt = ymax - ymin
	lzPt = zmax - zmin
	lMinPt = min(lxPt, lyPt, lzPt)

	# Deduced grid properties:
	spac = lMinPt / nMin
	dx = dy = dz = spac
	Nx, Ny, Nz = ceil(lxPt / spac), ceil(lyPt / spac), ceil(lzPt / spac)  # these are gridpoints strictly fitting the cloud
	nGpToUse = Vector3i(Nx, Ny, Nz) + Vector3i(2, 2, 2)  # these will be the points of the actual RegularGrid, extended by one gp on each side:
	gridMin = Vector3(xmin, ymin, zmin) - spac * Vector3(1, 1, 1)
	gridMax = Vector3(xmax, ymax, zmax) + spac * Vector3(1, 1, 1)

	# Initialization of numpy arrays (to phi=inf, for instance):
	phi = inf * numpy.ones(nGpToUse)  # every gp is initially considered to be outside, by default
	inside = numpy.zeros(nGpToUse)

	# Detecting and modifying surface gridpoints/cells to phi=0:
	for l in range(len(ptsX)):
		# In which (i,j,k) cell is the point of index l:
		i = int((ptsX[l] - gridMin[0]) / dx)
		j = int((ptsY[l] - gridMin[1]) / dy)
		k = int((ptsZ[l] - gridMin[2]) / dz)
		# assigning 0 as a distance value to that cell: (TODO: could be improved)
		phi[i, j, k] = 0
		# TODO: we may work here several times on the same cell

	# Looking at the center, considered to be inside, to initialize the algorithm with some phi = -inf:
	center = (gridMin + gridMax) / 2.
	i = int((center - gridMin)[0] / dx)
	j = int((center - gridMin)[1] / dy)
	k = int((center - gridMin)[2] / dz)
	phi[i, j, k] = -inf  # assigning phi=-inf to center which is assumed to be inside
	inside[i, j, k] = 1

	# Loop to determine all inner points, starting from center and looping over neighbours:
	startCell = [(i, j, k)]
	while startCell != []:  # startCell (with mutable values) will be modified below
		i = startCell[0][0]
		j = startCell[0][1]
		k = startCell[0][2]
		if i == 0 or i == nGpToUse[0] - 1 or j == 0 or j == nGpToUse[1] - 1 or k == 0 or k == nGpToUse[2] - 1:
			# setInside has leaked until grid boundaries without seeing the external surface
			raise BaseException(
			        'Impossible to proceed: grid is too fine (your choice of nMin =', nMin, 'would give ~', nMin**2,
			        'surface gp) wrt input cloud (which has only', len(ptsX),
			        'surface points) or surface is not closed. Please decrease nMin below ~',
			        len(ptsX)**0.5
			)
			break
		phi, inside = setInside(i, j, k, phi, inside, startCell)

	# Forming the final grid:
	gridBis = RegularGrid(min=gridMin, spacing=spac, nGP=nGpToUse)
	if ((gridBis.max() - gridMax) / spac).maxAbsCoeff() > 1 + 1.e-7:
		raise BaseException('Error', (gridBis.max() - gridMax) / spac)
	return (phi, gridBis)


def setInside(i, j, k, phi, inside, startCell):
	'''Developer function for phiIniFromPts. Propagates the "inside" nature, when appropriate
    :param i: x index of starting point
    :param j: y index of starting point
    :param k: z index of starting point
    :param phi: distance grid
    :param inside: grid with same size than phi, full of 0 initially, some values will be set to 1 when we go through inner gridpoints
    :param startCell: a liste of indices corresponding to those gridpoints whose neighbours are still to check
    '''
	startCell.remove((i, j, k))
	for dx in [-1, 1]:
		inside[i + dx, j, k] = 1
		if not isfinite(phi[i + dx, j, k]) and phi[i + dx, j, k] > 0:
			phi[i + dx, j, k] = -inf
			ngbrCell = (i + dx, j, k)
			startCell.append(ngbrCell)
	for dy in [-1, 1]:
		inside[i, j + dy, k] = 1
		if not isfinite(phi[i, j + dy, k]) and phi[i, j + dy, k] > 0:
			phi[i, j + dy, k] = -inf
			ngbrCell = (i, j + dy, k)
			startCell.append(ngbrCell)
	for dz in [-1, 1]:
		inside[i, j, k + dz] = 1
		if not isfinite(phi[i, j, k + dz]) and phi[i, j, k + dz] > 0:
			phi[i, j, k + dz] = -inf
			ngbrCell = (i, j, k + dz)
			startCell.append(ngbrCell)
	return phi, inside
