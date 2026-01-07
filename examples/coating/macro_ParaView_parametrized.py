#### import the simple module from the paraview
from paraview.simple import *
import numpy as np
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
# find source
my_source = FindSource('ProgrammableFilter1')
### extra function for paramatrizing input 

def input_with_default(prompt, default=None):
    user_input = input(f"{prompt} [{str(default)}]: ")
    if user_input == '':
        return default
    return user_input

######### SETTINGS

nPhi = input_with_default("Enter the phi resolution (integer)", default=8)
nTheta  = input_with_default("Enter the theta resolution (integer)", default=4)
scaleMin  = input_with_default("Enter minimum thickness in color scale (float)", default=0.)
scaleMax  = input_with_default("Enter maximul thickness in color scale (float)", default=1)

nPhi = int(nPhi)
nTheta  = int(nTheta)
scaleMin  = float(scaleMin)
scaleMax  = float(scaleMax)

#granice przedziałów theta policzone na podstawie H
dH = np.linspace(-1,1,nTheta+1)
thetaLimRad = np.arccos(dH) 
thetaLimDeg = np.degrees(thetaLimRad)


for thetaNo in range(nTheta):
    for phiNo in range(nPhi):
        currentName = 'theta_{:d}_phi_{:d}'.format(thetaNo,phiNo)
        ############
        # create a new 'Calculator' / actually to just copy data
        dataCopy = Calculator(registrationName=currentName, Input=my_source)
        # Properties modified on calculator1
        dataCopy.Function = ''
        # get active view
        renderView1 = GetActiveViewOrCreate('RenderView')
        # show data in view
        currentDisplay = Show(dataCopy, renderView1, 'GeometryRepresentation')
        # trace defaults for the display properties.
        currentDisplay.Representation = 'Surface'
        # hide data in view
        Hide(my_source, renderView1)
        # show color bar/color legend
        currentDisplay.SetScalarBarVisibility(renderView1, False)
        # update the view to ensure updated data information
        renderView1.Update()
        # get color transfer function/color map for 'radius'
        radiusLUT = GetColorTransferFunction('radius')
        # get opacity transfer function/opacity map for 'radius'
        radiusPWF = GetOpacityTransferFunction('radius')
        # get 2D transfer function for 'radius'
        radiusTF2D = GetTransferFunction2D('radius')
        RenameProxy(dataCopy, 'sources', currentName)
        # rename source object
        RenameSource(currentName, dataCopy)
        # change representation type
        currentDisplay.SetRepresentationType('3D Glyphs')
        # Properties modified on currentDisplay
        currentDisplay.Orient = 1
        # Properties modified on currentDisplay
        currentDisplay.OrientationMode = 'Quaternion'
        currentDisplay.SelectOrientationVectors = 'Quaternion'
        # Properties modified on currentDisplay
        currentDisplay.Scaling = 1
        # Properties modified on currentDisplay
        currentDisplay.ScaleFactor = 2.0
        currentDisplay.ScaleMode = 'Magnitude'
        currentDisplay.SelectScaleArray = 'radius'
        # Properties modified on currentDisplay
        currentDisplay.GlyphType = 'Sphere'

        ##### set some sphere settings - NOTE, THAT YADE HAS THE NOTATION DIFFERENT FROM PARAVIEW (PHI AND THETA ARE SWITCHED)
        currentDisplay.GlyphType.StartTheta = phiNo*360/nPhi
        currentDisplay.GlyphType.EndTheta = (phiNo+1)*360/nPhi
        currentDisplay.GlyphType.StartPhi = thetaLimDeg[thetaNo+1]
        currentDisplay.GlyphType.EndPhi = thetaLimDeg[thetaNo]
        # set scalar coloring
        ColorBy(currentDisplay, ('POINTS', currentName))

        # Hide the scalar bar for this color map if no visible data is colored by it.
        HideScalarBarIfNotNeeded(radiusLUT, renderView1)

        # rescale color and/or opacity maps used to include current data range
        currentDisplay.RescaleTransferFunctionToDataRange(True, False)

        # show color bar/color legend
        currentDisplay.SetScalarBarVisibility(renderView1, False)

        # get color transfer function/color map for 
        LUT = GetColorTransferFunction(currentName)

        # get opacity transfer function/opacity map 
        PWF = GetOpacityTransferFunction(currentName)

        # get 2D transfer function
        TF2D = GetTransferFunction2D(currentName)
        
        ##
        # Rescale transfer function
        LUT.RescaleTransferFunction(scaleMin, scaleMax)

        # Rescale transfer function
        PWF.RescaleTransferFunction(scaleMin, scaleMax)
        
        # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
        LUT.ApplyPreset('just_blue', True)
        # Apply a preset using its name.
        PWF.ApplyPreset('just_blue', True)

        # get layout
        layout1 = GetLayout()
