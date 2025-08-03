# -*- encoding=utf-8 -*-
# Vasileios Angelidakis <v.angelidakis@qub.ac.uk>
# Visco-elastic contact law using MultiScGeom and MultiViscElPhys, for multiple contacts per pair of particles

# The normal and shear stiffness and coefficient of restitution values can be defined either at the material level, or directly in the Ip2

#---------------------------------------------------------------------
# Material
O.materials.append(ViscElMat(young=-1, kn=4000, ks=4000, en=0.6, et=0.6, poisson=0.4, density=1000, frictionAngle=atan(0.5), label='particles')) 

#---------------------------------------------------------------------
# Free body
O.bodies.append(levelSetBody('sphere', (0,0,0.03), 0.01, spacing=0.0005, material='particles', nodesPath=1,nSurfNodes=83))
O.bodies[-1].shape.color=[0.2,0.2,0.8]
#O.bodies[-1].state.ori=Quaternion((1,0,0),pi/2) # Change the particle orientation to make a different number of nodes interact

#---------------------------------------------------------------------
# Fixed body
O.bodies.append(levelSetBody('sphere', (0,0,0), 0.01, spacing=0.0005, material='particles' ,nodesPath=1, nSurfNodes=83))
O.bodies[-1].shape.color=[0.8,0.2,0.2]
O.bodies[-1].state.blockedDOFs='xyzXYZ'
#O.bodies[-1].state.ori=Quaternion((1,0,0),pi/2)

#---------------------------------------------------------------------
# Engines
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_LevelSet_Aabb()], verletDist=0.0, label="collider"),
	InteractionLoop(
		[Ig2_LevelSet_LevelSet_MultiScGeom(label='ig2')],
		[Ip2_ViscElMat_ViscElMat_MultiViscElPhys(label='ip2')], # kn=4000, ks=4000, en=0.6, et=0.6 can be defined directly here too, instead of the material
		[Law2_MultiScGeom_MultiViscElPhys_Basic(label='law')],
		label='ILoop'),
	NewtonIntegrator(damping=0.0, gravity=(0,0,0), label='newton'), # no gravity
	]

O.dt=1e-6

#---------------------------------------------------------------------
# Set initial velocity before collision
b = O.bodies[0]
v0 = -0.4
b.state.vel = [0,0,v0]

#---------------------------------------------------------------------
# Plot results
from yade import plot
def addPlotData():
	plot.addData(	iteration	=	O.iter, 
					time 		=	O.time, 
					velocity	=	b.state.vel.norm(), # magnitude of velocity
#					velZ 		=	abs(b.state.vel[2]  # z component of velocity
				)

	# Pause simulation after the free particle rebounces
	if O.bodies[0].state.pos[2] > 0.03:
		O.pause()

plot.plots = {'iteration': ('velocity')} #,'velZ'
plot.plot()
O.engines=O.engines+[PyRunner(command='addPlotData()',iterPeriod=1, dead=False, label='addData')]

#---------------------------------------------------------------------
# Visualisation
from yade import qt
Gl1_LevelSet.surfNodes = True	# Render surface nodes
Gl1_LevelSet.wire = True		# Render wireframe
Gl1_LevelSet.surface = True 	# Render particle surface

v = qt.View()					# Activate qt view

v.eyePosition=Vector3(0,-0.075,0.015)
v.upVector=Vector3(0,0,1)
v.viewDir=Vector3(0,1,0)

#---------------------------------------------------------------------
# Run enough timesteps for a contact to happen, and compare the requested coefficient of restitution with the observed one.
O.run(50000,True)
v1 = O.bodies[0].state.vel[2]

print('\nThe requested (input) normal coefficient of restitution was', O.materials[0].en)
print('The observed (output) normal coefficient of restitution is', abs(v1/v0))
print('The relative error between the two is', abs(abs(v1/v0)-O.materials[0].en)/O.materials[0].en,'\n')

