These scripts aim to illustrate and/or facilitate the use of the [LevelSet](https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.LevelSet) class for shape description.

Please note in particular the following simulation examples:

- `discharge.py` for a LevelSet-based simulation of a discharge of 1000 superquadric particles, as discussed in Section 6 of [[Duriez2021b](https://www.sciencedirect.com/science/article/pii/S0098300421002247)]
- `levelSetBody.py` for a particle-scale (1 or 2 bodies) simulation in order to illustrate the various syntaxes available for defining a LevelSet-shaped body in YADE and the expected `O.engines`
- `lsNodeGeom.py` for illustrating the use of [Ig2_LevelSet_LevelSet_LSnodeGeom](https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.Ig2_LevelSet_LevelSet_LSnodeGeom) for faster simulations in the case of a particle-scale example of a bouncing LevelSet sphere on a LevelSet floor
- `multiScGeom.py` for a particle-scale example of the [MultiScGeom](https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.MultiScGeom) description of a contact between non-convex LevelSet shapes (with several contact points)
- `rendering.py` for a particle-scale example of rendering options

The other files provide some utilities functions, each in a separate file, a number of those being dedicated to create LevelSet-shaped bodies from .stl shape data:

- `mkStlInertial.py` for transforming a .stl surface mesh into an inertial configuration as this is expected by `phiFromStl` or `phiIniFromPts`
- `phiFromStl.py` for obtaining a distance field from a .stl, thanks in particular to the external [trimesh](https://pypi.org/project/trimesh/) Python module. See also `phiIniFromPts`
- `phiIniFromPts.py` for the same purpose than `phiFromStl`, with a required extra step of a Fast Marching Method to perform by the user, more zeros in the distance field but lighter computational costs and no external dependancy
- `pvVisu.py` for a Python-automated vizualization of LevelSet-shaped YADE bodies in Paraview (after using `VTKRecorder` in the simulation, see `levelSetBody.py` for an example)
- `ramFP.py` for measuring current RAM usage on your machine, as this could be a limiting factor when using LevelSet
- `spaceRot.py` for defining "random" orientations in 3D space