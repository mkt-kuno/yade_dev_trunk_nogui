import numpy as np
import matplotlib.pyplot as plt

# Read the data from the file
data = np.loadtxt('TimeData.txt')

# Filter rows where the first column (iterations) is larger than 200,000
filtered_data = data[data[:, 0] > 200000]

# Extract columns from filtered data
Iteration = filtered_data[:, 0]  # First column
Height = filtered_data[:, 1] # Second column
Phi = filtered_data[:, 2] # Third column

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(Iteration, Height/Height[0])

# Adding labels and title
plt.xlabel('Iteration')
plt.ylabel('$H/H_0$')
plt.title('Height as a Function of Iterations')
plt.savefig('HeightEvolution.png')
plt.show()


# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(Iteration, Phi, color='orange')

# Adding labels and title
plt.xlabel('Iteration')
plt.ylabel('Volume fraction')
plt.title('Phi as a Function of Iterations')
plt.savefig('phiEvolution.png')
plt.show()
