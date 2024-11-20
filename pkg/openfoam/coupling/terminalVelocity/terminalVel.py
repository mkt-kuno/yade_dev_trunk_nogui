import matplotlib.pyplot as plt
import numpy as np
def read_data(filename):
    times = []
    velocities = []
    positions = []

    with open(filename, 'r') as file:
        # Skip header line if there is one
        header = file.readline()

        for line in file:
            # Split the line into values
            values = line.strip().split()
            if len(values) == 3:
                time, velocity, position = map(float, values)
                times.append(time)
                velocities.append(velocity)
                positions.append(position)

    return times, velocities, positions

def read_forces(filename):
    Fx = []
    Fy = []
    Fz = []

    with open(filename, 'r') as file:
        for line in file:
            values = line.strip().split()
            if len(values) == 3:
                fx, fy, fz = map(float, values)
                Fx.append(fx)
                Fy.append(fy)
                Fz.append(fz)
    return Fx, Fy, Fz


def smooth(data, window_size=50):
    """Apply a moving average to smooth data."""
    return np.convolve(data, np.ones(window_size)/window_size, mode='same')


################ Analytical solution ############
densityF=1000
nu=0.0001
mu=densityF*nu
densityParticle=1500
Radius=0.0005
g=-9.81
Vol=4*np.pi*Radius**3/3

weightValues=[densityParticle*Vol*g,densityParticle*Vol*g]
ArchValues=[densityF*Vol*-g,densityF*Vol*-g]


TermVel_analytical=(2*Radius*Radius*(densityParticle-densityF)*g)/(9*mu)
Reyonolds=abs(densityF*2*Radius*TermVel_analytical)/mu
print ("Reynolds number: ",Reyonolds)
print ("Terminal velocity: ",TermVel_analytical)


################ Plots ############

filename = 'data.txt'
times, velocities, positions = read_data(filename)

filename_forces = 'archimedes_force.txt'
Fx_arch, Fy_arch, Fz_arch = read_forces(filename_forces)
filename_forces = 'hydro_drag_force.txt'
Fx_drag, Fy_drag, Fz_drag = read_forces(filename_forces)

# Smooth forces
Fy_arch_smooth = smooth(Fy_arch)
Fy_drag_smooth = smooth(Fy_drag)

# Calculate Total Force in Y direction
Fy_Total = [fy_drag + fy_arch + densityParticle * Vol * g
            for fy_drag, fy_arch in zip(Fy_drag_smooth, Fy_arch_smooth)]

timeStep=1e-6
maxT= (len(Fy_drag)*timeStep)
timesF=np.arange(0,maxT,timeStep)


timeAnalytical=[times[0],times[-1]]
velAnalytical=[TermVel_analytical,TermVel_analytical]


fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 12))
# Plot Time vs Velocity
axs[0].plot(times, velocities, marker='None', linestyle='--', color='b')
axs[0].plot(timeAnalytical, velAnalytical, marker='None', linestyle='-', color='k')
# axs[0].set_xlabel('Time')
axs[0].set_ylabel('Velocity')

# Plot Time vs Position
axs[1].plot(times, positions, marker='None', linestyle='--', color='b')
axs[1].set_ylabel('Position')


# Plot Forces
axs[2].plot(timesF, Fy_arch_smooth, label='Archemedes force', color='b')
axs[2].plot(timesF, Fy_drag_smooth, label='Drag force', color='g')
axs[2].plot(timeAnalytical, weightValues, label='Weight', color='r')
# axs[2].plot(timeAnalytical, ArchValues, label='Analytical archemedes', color='navy')
axs[2].plot(timesF, Fy_Total, label='Total force', color='k')

axs[2].set_xlabel('Time [s]')
axs[2].set_ylabel('Force [N]')
axs[2].legend()

# Adjust layout
plt.tight_layout()

# Show plots
plt.savefig("terminalVel.png")
plt.show()
