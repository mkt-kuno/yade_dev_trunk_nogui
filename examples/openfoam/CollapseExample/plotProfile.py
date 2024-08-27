import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

# Read the points from the file
def read_points(file_path):
    points = []
    with open(file_path, 'r') as file:
        for line in file:
            points.append(list(map(float, line.split())))
    return np.array(points)

# Find the upper contour based on the y component
def find_upper_contour_y(points,xStep):
    # Use only x and y components
    points_2d = points[:, :2]
    # Sort points by x-coordinate
    sorted_points = points_2d[np.argsort(points_2d[:, 0])]
    # Find the upper contour
    upper_contour = []

    lenX=np.linspace(0, sorted_points[-1][0], int(sorted_points[-1][0]/xStep))
    for i in range(len(lenX)):
       
        max_y = 0
        max_x = lenX[i]

        for point in sorted_points:

            if (point[1] > max_y) and (abs(point[0] -lenX[i]) < xStep):
                max_y = point[1]
                max_x = point[0]
        upper_contour.append((max_x,max_y))

    return np.array(upper_contour)

# Main function

directory = 'data'

# Get a list of all files in the directory
files = os.listdir(directory)

# Regular expression to match the numerical part of the filenames
pattern = re.compile(r'profile_(\d+)\.txt')

# Extract numerical parts from the filenames
numbers = [int(pattern.search(f).group(1)) for f in files if pattern.search(f)]

# Find the minimum and maximum values
Ini = min(numbers)
Final = max(numbers)
timeStep=1e-5
Nsteps=10000
time = int(timeStep * Final * 1e8) / 1e8  # Multiply and divide to maintain precision
Radius=3.84e-3/2.*2
xStep=Radius*1
Lini=0.04


value = round(time / timeStep)



file_path = directory+'/profile_%i.txt'%value
points = read_points(file_path)
upper_contour_2d = find_upper_contour_y(points,xStep)
#plot_upper_contour_y(points, upper_contour_2d)


X_exp_a2=[]
Y_exp_a2=[]
with open("expSun_a2.txt", "r") as file:
    next(file)  # Skip the first line containing column labels
    for line in file:
        columns = line.split()
        X_exp_a2.append(float(columns[0]))
        Y_exp_a2.append(float(columns[1]))


        
        

# Plot all points
plt.scatter(points[:, 0], points[:, 1], c='blue', alpha=0.1, marker='o')
# Plot the upper contour
plt.plot(upper_contour_2d[:, 0], upper_contour_2d[:, 1], c='red',label='DEM-OF')
plt.plot(X_exp_a2, Y_exp_a2, c='k',marker='s',linestyle='none',label='Exp. Sun (2020)')

plt.xlabel('X [m]')
plt.ylabel('Y [m]')
plt.xlim([0,0.2])
plt.legend()
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
#plt.title('Upper Contour of Points (Based on Y)')
plt.savefig('ProfilePerm_%s.png'%directory)
plt.show()
