# -*- encoding=utf-8 -*-
# Vasileios Angelidakis <v.angelidakis"qub.ac.uk> Jerome Duriez <jerome.duriez@inrae.fr>
# Rendering options of LevelSet particles

from yade import qt

O.bodies.append(levelSetBody('sphere',(3,0,0),1))
O.bodies[-1].shape.color=[0.2,0.9,0.3]

Gl1_LevelSet.surfNodes = True	# Render surface nodes
Gl1_LevelSet.wire = True		# Render wireframe
Gl1_LevelSet.surface = True 	# Render particle surface

v = qt.View()					# Activate qt view
