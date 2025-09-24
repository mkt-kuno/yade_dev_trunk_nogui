"""
This file implements a number of clump algorithms described by Fathipour-Azar & Duriez (2025). 

Reference:
Fathipour-Azar, H., Duriez, J. A Comparative Study of Existing and New Sphere Clump Generation Algorithms for Modeling Arbitrary Shaped Particles. Arch Computat Methods Eng (2025). https://doi.org/10.1007/s11831-025-10256-1
"""

try:
    import trimesh
    import numpy as np
    import scipy.ndimage as nd
    import csv
    from sklearn.cluster import KMeans
except:
    raise BaseException('One dependency is not required, please check the doc (readme.md) and your operating system')

##################################################
#### Auxiliary functions for playing with .stl ###
##################################################

def binarize(filename, voxel_size):
    '''Returns a 3D numpy array serving as a binary representation (1 = inside, 0 = outside) of the inner volume of the surface described by *filename*
    :param str filename: path to the file surface at hand
    :param float voxel_size: which voxel size to use in the binary representation
    '''
        
    # Load the file using the appropriate method based on the file type
    
    # If the input is an STL file, it will be loaded using the trimesh.load method:
    if filename.endswith('.stl'):
        mesh = trimesh.load(filename, file_type='stl')
        
    # If the input is a PLY file, it will be loaded using the trimesh.load_ply method.    
    elif filename.endswith('.ply'):
        mesh = trimesh.load_ply(filename, file_type='ply')
        
    # If the input is a point cloud, you can create a mesh from it using the `points` method:
    elif filename.endswith('.xyz'):
        points = trimesh.points.load(filename) #"path/to/file.xyz"
        mesh = trimesh.Trimesh(vertices=points) # , process=True)
        
    # If the input is a NumPy .npy file, it will be loaded using np.load and then passed to the trimesh.Trimesh constructor:
    # Load 3D binary image
    elif filename.endswith('.npy'):
        points = np.load(filename, file_type='npy')
        mesh = trimesh.Trimesh(vertices=points)
        #mesh = trimesh.points.PointCloud(points)
    # Add other file types as needed
    else:
        raise BaseException('Only .stl, .ply, .xyz and .npy are presently supported for *filename*')
    
    # After voxelizing the above mesh data (just the surface !) for a given voxel size through .voxelized (ie, turning it into a trimesh.VoxelGrid instance, whose .matrix is an expected 3D numpy.array made of True / False)
    # The fill method is then called on the voxelized_mesh object to fill in the interior of the mesh with voxels.
    voxelized_mesh = mesh.voxelized(pitch=voxel_size).fill()
    
    # Set the inside voxels (where True holds) to 1 and the outside voxels to 0
    voxels = np.where(voxelized_mesh.matrix, 1, 0)
    
    # Save the voxels as a mesh
    #voxels.export('voxels.ply')

    return voxels


def get_stl_stats(stl_path):
    '''Returns volume, surface and inertia for *stl_path* (for any kind of user-purpose)'''
    # Load the STL file into a trimesh mesh object
    mesh = trimesh.load(stl_path)

    # Calculate the volume of the mesh
    volume = mesh.volume

    # Calculate the surface area of the mesh
    surface_area = mesh.area

    # Calculate the inertia tensor of the mesh
    inertia = mesh.moment_inertia

    return volume, surface_area, inertia

### End of .stl-oriented functions ###

#########################################
#### Actual clump algorithm functions ###
#########################################

def distance_transform(binary_image, voxel_size, threshold, out_csv = './Data/Detailed_medSurf_spheres.csv'):
    '''Applies the "Detailed aproach" of the manuscript § 2.1.1 to *binary_image*, ie returns the corresponding medial surface distance-transform values in a 3D numpy array, after creating a .csv with the information of a detailed clump version, depending on *out_csv*
    :param 3D numpy array binary_image: the binary image (made of 0 and 1) serving as an input
    :param float voxel_size: the underlying voxel size to *binary_image* to have correct real-world units in *out_csv*
    :param float threshold: enables to disregard (i.e. consider as void) from the returned object and *out_csv* the medial surface voxels that are closer to the boundary than *threshold* (in voxels unit) and that would correspond (during clump sphere creation) to insignificant spheres, with radii smaller than *threshold*. A 0 value makes it ineffective
    :param str out_csv: the path where to store the detailed clump .csv in x,y,z,r format. Clump data obey real world units provided *voxel_size* is consistent to the origins of *binary_image* but its Aabb will always start at the space origin. *out_csv* may include slashes but corresponding (sub-)folders have to exist already. No clump .csv is stored if empty.
    :returns: a 3D numpy array of medial surface distance transform (dt) values in voxels unit for the binary image (NB: in the process, all "outside" voxels are assigned a zero dt value. See also https://stackoverflow.com/a/44770662)
    '''
    
    # Add zero band around binary image using pad function
    margin = 1
    v_margin = np.pad(binary_image, [(margin, margin), (margin, margin), (margin, margin)], mode='constant')

    # Initialize medial_surface array with the shape of binary_image
    medial_surface = np.zeros(shape=binary_image.shape, dtype='i4')
    
    # Calculate the distance map using the Euclidian distance transform
    sedt_margin = (nd.distance_transform_edt(v_margin) ).astype('i4')
    
    # Create a 3D kernel for the non-maximum suppression algorithm
    kernel = np.array([
        [[0, 0, 0], [0, 1, 0], [0, 0, 0]],
        [[0, 1, 0], [1, 1, 1], [0, 1, 0]],
        [[0, 0, 0], [0, 1, 0], [0, 0, 0]]
    ])

    # Perform non-maximum suppression on the distance map
    medial_surface = nd.maximum_filter(sedt_margin, footprint=kernel)
    medial_surface[medial_surface != sedt_margin] = 0

    # Remove zero band around the binary image
    medial_surface = medial_surface[margin*2:-margin*2, margin*2:-margin*2, margin*2:-margin*2]
    
    # check threshold 
    if threshold < medial_surface.min() or threshold > medial_surface.max():
        raise ValueError("Given threshold is " + str(threshold) + " but should be between the min = " + str(medial_surface.min()) + " and the max = " + str(medial_surface.max()) + " of the distance transform values for the input image")
    if threshold > 0: # no need to apply the below line if threshold = 0, that would leave medial_surface strictly unchanged
        medial_surface[medial_surface < threshold] = 0
    
    if out_csv: # condition is False iff out_csv is empty
        # Open the CSV file for writing
        with open(out_csv, 'w') as csv_file:
            # Write the header line
            csv_file.write('#x_center\ty_center\tz_center\tradius\n')
            # Iterate over the voxels in the grid
            for i in range(medial_surface.shape[0]):
                for j in range(medial_surface.shape[1]):
                    for k in range(medial_surface.shape[2]):
                        if medial_surface[i,j,k] > 0:
                            csv_file.write(f'{i*voxel_size}\t{j*voxel_size}\t{k*voxel_size}\t{medial_surface[i,j,k]*voxel_size}\n')
    
    # Return the distance transform values (as a 3D numpy array)
    return medial_surface

def volCov(binary_image, voxel_size, threshold, coverage = 1, out_csv = './Data/Selected_medSurf_spheres.csv', verbose = 0):
    """
    Corresponding to the Greedy volume coverage § 2.1.2 of the manuscript, the function proposes a clump description of *binary_image*, stored in *out_csv*, by working out volume coverage considerations on the medial surface distance transform of the *binary_image*
    :param 3D numpy array binary_image: the binary image (made of 0 and 1) serving as an input
    :param float voxel_size: the underlying voxel size to *binary_image*, to have correct real-world units in *out_csv*
    :param float threshold: to set up (in voxels unit) a minimum radius value for the sphere members, after passing it to the corresponding argument of distance_transform. A 0 value avoids all filtering
    :param float coverage: a ]0; 1] proportion of voxels that will be considered for creating spheres. Serving as a proxy for volume precision
    :param str out_csv: path where to store the detailed clump .csv in x,y,z,r format. Clump data obey real world units provided *voxel_size* is consistent to the origins of *binary_image* but its Aabb will always start at the space origin. *out_csv* may include slashes but corresponding (sub-)folders have to exist already.
    :param bool verbose: whether to print (if True) information messages during the iterative process
    :returns: nothing
    """
    import csv
    import numpy as np

    # Compute medial surface dt for binary image
    medial_surface = distance_transform( binary_image, voxel_size, threshold, '')

    # Find the non-zero voxels of binary_image
    non_zero_voxels_bi = np.transpose(np.nonzero(binary_image))

    # Find the non-zero voxels of distance_map
    non_zero_voxels_dmap = np.transpose(np.nonzero(medial_surface))

    # Initialize a variable to store the sum of voxels
    coverage_sum = np.zeros(len(non_zero_voxels_dmap))

    # get the number of initial non zero voxels in binary image
    initial_non_zero_voxels = len(non_zero_voxels_bi)

    # open a csv file to save the selected distance map voxels in each iteration
    with open(out_csv, 'w') as csvfile:
        writer = csv.writer(csvfile,delimiter = '\t')
        # write header
        writer.writerow(['#x_center','y_center','z_center','radius', 'volume'])

        max_coverage = 1 # just initializing a below-loop-variable to any strictly positive value

        itCptr = 0
        while len(non_zero_voxels_bi) > initial_non_zero_voxels*(1-coverage) and max_coverage > 0:
            if verbose:
                print('After',itCptr,'iterations, we still have',(100.*len(non_zero_voxels_bi))/initial_non_zero_voxels,'percent of given solid voxels to handle, starting a new one')
            coverage_sum = np.zeros(len(non_zero_voxels_dmap))

            # Specific case where the distance matrix is too large for the memory, If it is, it splits the non-zero voxels of the binary image into smaller chunks, and for each chunk, it calculates the euclidean distance between it and the non-zero voxels of the distance map. Then it adds the result of each chunk together.
            if len(non_zero_voxels_bi) * len(non_zero_voxels_dmap) > 200e6:
                n = max(1, int(len(non_zero_voxels_bi) * len(non_zero_voxels_dmap) / 50e6))
                print('Splitting the data into',n,'chunks, to handle one after another for computational reasons')
                non_zero_voxels_bi_s = np.array_split(non_zero_voxels_bi, n)

                for j in range(n):
                    non_zero_voxels_bi_tmp = non_zero_voxels_bi_s[j]

                    # Calculate euclidean distance between x,y,z coordination of voxels
                    euclidean_distance = np.linalg.norm(non_zero_voxels_dmap[:, None] - non_zero_voxels_bi_tmp, axis=-1)

                    # check if euclidean distance is smaller than the value of that voxel in distance map
                    mask = euclidean_distance < medial_surface[non_zero_voxels_dmap[:, 0],non_zero_voxels_dmap[:, 1],non_zero_voxels_dmap[:, 2]][:,None]
                    coverage_sum += np.sum(mask, axis=1)

            else:
                # Calculate euclidean distance between x,y,z coordination of voxels
                euclidean_distance = np.linalg.norm(non_zero_voxels_dmap[:, None] - non_zero_voxels_bi, axis=-1)

                # check if euclidean distance is smaller than the value of that voxel in distance map
                mask = euclidean_distance < medial_surface[non_zero_voxels_dmap[:, 0],non_zero_voxels_dmap[:, 1],non_zero_voxels_dmap[:, 2]][:,None]
                coverage_sum = np.sum(mask, axis=1)

            # find the max coverage_sum
            max_coverage = max(coverage_sum)

            # find the index of max coverage_sum
            max_index = np.argmax(coverage_sum)

            # find the x,y,z and value of the distance map voxel that covers the most
            x, y, z = non_zero_voxels_dmap[max_index]
            value = medial_surface[x, y, z]

            # save x,y,z and value of the selected voxel in distance map to a csv file
            writer.writerow([x*voxel_size, y*voxel_size, z*voxel_size, value*voxel_size, max_coverage])

            # Set binary image voxels that are covered by the max distance map voxel to 0
            # Generator : the loop only iterates over the voxels that are currently non-zero, rather than iterating over all the voxels in the array.
            #for bi_voxel in (i for i in non_zero_voxels_bi):
            #    euclidean_distance = np.linalg.norm(non_zero_voxels_dmap[max_index] - bi_voxel)
            #    if euclidean_distance < medial_surface[x, y, z]:
            #        binary_image[bi_voxel[0], bi_voxel[1], bi_voxel[2]] = 0
            euclidean_distance = np.linalg.norm(non_zero_voxels_dmap[max_index] - non_zero_voxels_bi, axis=-1)
            mask = euclidean_distance < medial_surface[x, y, z]
            binary_image[non_zero_voxels_bi[mask, 0], non_zero_voxels_bi[mask, 1], non_zero_voxels_bi[mask, 2]] = 0

            # update non_zero_voxels_bi
            non_zero_voxels_bi = np.transpose(np.nonzero(binary_image))

            itCptr += 1
            if verbose:
                print('At the end of iteration',itCptr,', we still have', (binary_image > 0).sum(),'voxels of solid to handle (when looking at a .sum(), or',len(non_zero_voxels_bi),'when looking at non_zero_voxels_bi)')
                print('Which will be compared with',initial_non_zero_voxels*(1-coverage),'to decide whether we go for another round. See also max_coverage =',max_coverage,'which has to be strictly positive')
            # Terminate the loop when the remaining non-zero voxels in the binary image is less than x of the initial non-zero voxels or when we stopped gaining anything in volume coverage.

def clustKMeans(binary_image, voxel_size, n_clusters, smallestToKeep = 0, out_csv = './Data/Clustered_spheres.csv'):
    """
    Corresponding to the K‐means Clustering Algorithm § 2.1.3 of the manuscript, the code uses the KMeans method from the scikit-learn library to cluster the clump data from an initial detailed clump model of *binary_image*, into just *n_clusters* clump components
    :param 3D numpy array: the binary image (made of 0 and 1) serving as an input
    :param float voxel_size: the underlying voxel size to *binary_image* to have correct real-world units in *out_csv*
    :param int n_clusters: how many clusters you wish to obtain through the inner KMeans approach, ie, how many sphere components in the clump (modulo *smallestToKeep* considerations)
    :param float smallestToKeep: use a strictly positive value to disregard from the final clump description, initial clump members with a radius smaller than *smallestToKeep* times absolute minimum of radius. Use any negative value (e.g., 0) to keep all clump members. On the other hand, if some filtering of the smallest spheres is desired for computational efficiency, 2 had been observed to give satisfactory results
    :param str out_csv: the path where to store the final clump .csv in x,y,z,r format. Clump data obey real world units provided *voxel_size* is consistent to the origins of *binary_image* but its Aabb will always start at the space origin. *out_csv* may include slashes but corresponding (sub-)folders have to exist already.
    """

    import csv
    from sklearn.cluster import KMeans

    distance_transform(binary_image, voxel_size, 0, '/tmp/Detailed.csv')

    # Load the data from the above CSV file in x y z r format and in real world units
    data = []
    with open('/tmp/Detailed.csv', 'r') as f:
        reader = csv.reader(f, delimiter = '\t')
        next(reader) # skip header row
        for row in reader:
            data.append([float(x) for x in row])

    # Cluster the data
    kmeans = KMeans(n_clusters=n_clusters, random_state=1001) # doc at https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html
    kmeans.fit(data)
    centers = kmeans.cluster_centers_ # x y z r

    # Define a threshold for the clump members radii
    threshold = smallestToKeep * min([center[3] for center in centers])

    # Save the clusters to a new CSV file
    written_centers = 0
    with open(out_csv, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['#x_center', 'y_center', 'z_center', 'radius'])
        for center in centers:
            if center[3] < threshold:
                continue
            writer.writerow(center)
            written_centers += 1
    if threshold > 0:
        print(f'Clustering for n_clusters = {n_clusters} done but {written_centers} clump members considering a threshold radius = {threshold} (because of given smallestToKeep = {smallestToKeep}) were saved in .csv.')
