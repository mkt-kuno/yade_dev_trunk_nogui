import numpy as np

def sph_to_cart(theta, phi):
    """
    Convert spherical coordinates to carthesian.
    """
    x = np.cos(theta) * np.cos(phi)
    y = np.cos(theta) * np.sin(phi)
    z = np.sin(theta)
    return Vector3(float(x), float(y), float(z))

# function preparing triangle pattern
def triangulateSegmentedSphere(nTheta=4, nPhi=8, quality = 1):
    """
    Triangulates sphere and prepares a mapper to append field value to each field. 
    If quality > 0, eqch field is divided in more triangles (by multiplying phi and theta resolution by (1+quality))
    Returns
    - points
    - triangles
    - cell mapper (list of segment positions related to the triangles)
    """
    assert nTheta > 1 and nPhi > 2
    
    assert type(quality) == int and quality >=0 

    # multiplicator to increase number of fields (traingles)
    multiplicator = (1 + quality)# multiply the number of triangles by this factor
    
    nTheta *= multiplicator
    nPhi *= multiplicator
    
    # theta angles (equal heights / areas)
    theta = [np.arcsin(-1 + 2*i/nTheta) for i in range(nTheta+1)]
    
    # phi angles (equal angles)
    phi = [2*np.pi*j/nPhi for j in range(nPhi+1)]
    
    points = []
    index_map = {}
    
    # points for fields in the middle (divided by meridians)
    for i in range(1, nTheta):
        for j in range(nPhi):
            p = sph_to_cart(theta[i], phi[j])
            index_map[(i, j)] = len(points)
            points.append(p)
    
    # points for fields near poles
    south_pole_idx = len(points); points.append(sph_to_cart(-np.pi/2, 0.0))
    north_pole_idx = len(points); points.append(sph_to_cart(+np.pi/2, 0.0))
    
    triangles = []
    cell_scalars_mapper = []
    
    def segmentNo(thetaNo, phiNo, nPhi):
        return (thetaNo // multiplicator) * (nPhi // multiplicator) + phiNo // multiplicator
    
    # north pole (only triangle fields)
    for phiNo in range(nPhi):
        a = south_pole_idx
        b = index_map[(1, phiNo)]
        c = index_map[(1, (phiNo+1)%nPhi)]
        triangles.append((a,b,c))
        #cell_scalars.append(float(fieldValues[segmentNo(0,phiNo,nPhi)]))
        cell_scalars_mapper.append(segmentNo(0,phiNo,nPhi))
        
    # meridian stripes // each field divided into two triangles
    for thetaNo in range(1, nTheta-1):
        for phiNo in range(nPhi):
            v00 = index_map[(thetaNo, phiNo)]
            v01 = index_map[(thetaNo, (phiNo+1)%nPhi)]
            v10 = index_map[(thetaNo+1, phiNo)]
            v11 = index_map[(thetaNo+1, (phiNo+1)%nPhi)]
            # dicision into triangles
            triangles.append((v00, v10, v11))
            #cell_scalars.append(float(fieldValues[segmentNo(thetaNo,phiNo,nPhi)]))
            cell_scalars_mapper.append(segmentNo(thetaNo,phiNo,nPhi))
            triangles.append((v00, v11, v01))
            #cell_scalars.append(float(fieldValues[segmentNo(thetaNo,phiNo,nPhi)]))
            cell_scalars_mapper.append(segmentNo(thetaNo,phiNo,nPhi))
            
    # south pole (only triangle fields)
    for phiNo in range(nPhi):
        a = index_map[(nTheta-1, (phiNo+1)%nPhi)]
        b = index_map[(nTheta-1, phiNo)]
        c = north_pole_idx
        triangles.append((a,b,c))
        #cell_scalars.append(float(fieldValues[segmentNo(nTheta-1,phiNo,nPhi)]))
        cell_scalars_mapper.append(segmentNo(nTheta-1,phiNo,nPhi))
    
    return points, triangles, cell_scalars_mapper

def write_vtk(filename, points, triangles, cell_scalars, field_name = 'coating_thickness'):
    """Simple vtk writer."""
    with open(filename, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Triangulated sphere\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        
        # points
        f.write(f"POINTS {len(points)} float\n")
        for x,y,z in points:
            f.write(f"{x} {y} {z}\n")
        
        # triangles
        n_tris = len(triangles)
        f.write(f"POLYGONS {n_tris} {n_tris*4}\n")
        for a,b,c in triangles:
            f.write(f"3 {a} {b} {c}\n")
        
        # scalars
        f.write(f"CELL_DATA {n_tris}\n")
        f.write("SCALARS {:s} float 1\n".format(field_name))
        f.write("LOOKUP_TABLE default\n")
        for s in cell_scalars:
            f.write(f"{float(s):.6f}\n")
    

## actual function exporting spheres
class SegmentedSphereExporter:
    def __init__(self, filebase, nTheta=4, nPhi=8, quality=1):
        self.filebase = filebase
        self.nTheta = nTheta
        self.nPhi = nPhi
        self.quality = quality
        self.counter = 0  # 

    def __call__(self):
        bodyCount = 0
        # precompute triangles pattern
        pts, tris, scalars_mapper = triangulateSegmentedSphere(
            self.nTheta, self.nPhi, self.quality
        )
        # lists for points of all bodies
        all_pts = []
        all_tris = []
        all_scalars = []
        for b in O.bodies:
            if isinstance(b.shape, Sphere) and isinstance(b.state, SegmentedState):
                pos = b.state.pos
                ori = b.state.ori
                r = b.shape.radius

                points = [(ori * pt) * r + pos for pt in pts]

                shift = len(points) * bodyCount
                triangles = [(tri[0] + shift, tri[1] + shift, tri[2] + shift) for tri in tris]

                thickness = b.state.coatingThickness
                scalars = [float(thickness[segmentNo]) for segmentNo in scalars_mapper]

                all_pts += points
                all_tris += triangles
                all_scalars += scalars
                bodyCount += 1

        filename = f"{self.filebase}-segmented_spheres_{self.counter:05d}.vtk"
        write_vtk(filename, all_pts, all_tris, all_scalars, field_name='coating_thickness')
        self.counter += 1
    
            
            
