#### import the simple module from the paraview
from paraview.simple import *
import os
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

def createImages(path1,path2,outputPath):
    isoFaces_1vtk = LegacyVTKReader(FileNames=[path1])
    isoFaces_2vtk = LegacyVTKReader(FileNames=[path2])
    
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1084, 814]
    
    # show data in view
    #isoFaces_844vtkDisplay = Show(isoFaces_844vtk, isoFaces_1vtk)
    #isoFaces_844vtkDisplay = Show(isoFaces_844vtk, isoFaces_2vtk)
    
    #renderView1.ResetCamera()
    isoFaces_1vtkvtkDisplay = Show(isoFaces_1vtk, renderView1)
    isoFaces_2vtkvtkDisplay = Show(isoFaces_2vtk, renderView1)
    
    
    # trace defaults for the display properties.
    isoFaces_1vtkvtkDisplay.ColorArrayName = [None, '']
    isoFaces_2vtkvtkDisplay.ColorArrayName = [None, '']

    isoFaces_1vtkvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    isoFaces_1vtkvtkDisplay.SelectOrientationVectors = 'None'
    isoFaces_1vtkvtkDisplay.ScaleFactor = 0.042050239443778996
    isoFaces_1vtkvtkDisplay.SelectScaleArray = 'None'
    isoFaces_1vtkvtkDisplay.GlyphType = 'Arrow'
    isoFaces_1vtkvtkDisplay.GaussianRadius = 0.021025119721889498
    isoFaces_1vtkvtkDisplay.SetScaleArray = [None, '']
    isoFaces_1vtkvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    isoFaces_1vtkvtkDisplay.OpacityArray = [None, '']
    isoFaces_1vtkvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    
    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    isoFaces_1vtkvtkDisplay.OSPRayScaleFunction.Points = [0.0003351210034452379, 0.0, 0.5, 0.0, 0.9996629953384399, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    isoFaces_1vtkvtkDisplay.ScaleTransferFunction.Points = [0.0003351210034452379, 0.0, 0.5, 0.0, 0.9996629953384399, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    isoFaces_1vtkvtkDisplay.OpacityTransferFunction.Points = [0.0003351210034452379, 0.0, 0.5, 0.0, 0.9996629953384399, 1.0, 0.5, 0.0]
    
    
    # create a new 'Slice'
   
    # show data in view
   
    renderView1.ResetCamera()
    #change interaction mode for render view
    
    # current camera placement for renderView1
    # current camera placement for renderView1
    renderView1.CameraPosition = [0.40737003460526466, 0.5281376093626022, 2.8328371376901633]
    renderView1.CameraFocalPoint = [0.40737003460526466, 0.5281376093626022, 0.45497142150998116]
    renderView1.CameraParallelScale = 0.6154369340437763
        
    # destroy legacyVTKReader2
#    Delete(isoFaces_1vtk)
#    del isoFaces_1vtk
#    
#    # destroy legacyVTKReader1
#    Delete(isoFaces_2vtk)
#    del isoFaces_2vtk
    

#os.chdir('S:\\Verification\\VoFLibary\\release\\Advection\\vortexShearedDisc')
#testVar = raw_input()
#os.chdir('S:\\Verification\\VoFLibary\\release\\Advection\\vortexShearedDisc')

# create a new 'Legacy VTK Reader'
#isoFaces_1vtk = LegacyVTKReader(FileNames=['S:\\Verification\\VoFLibary\\release\\Advection\\vortexShearedDisc\\plicRDFN\\hex\\N128Co0.5\\isoFaces\\isoFaces_844.vtk'])
#isoFaces_2vtk = LegacyVTKReader(FileNames=['S:\\Verification\\VoFLibary\\release\\Advection\\vortexShearedDisc\\plicRDFN\\hex\\N128Co0.5\\isoFaces\\isoFaces_1683.vtk'])
#os.chdir('S:\\Verification\\VoFLibary\\release\\Advection\\deformationSphere')
os.chdir('S:/Verification/VoFLibary/release/Advection/deformationSphere/')
print os.listdir('./plicRDFN')
algorithm =  './plicRDFN'

isoFacesPath = [ 'S:/Verification\VoFLibary/release/Advection/deformationSphere/plicRDFN/hex/N64Co0.5/isoFaces/isoFaces_102.vtk','S:/Verification/VoFLibary/release/Advection/deformationSphere/plicRDFN/hex/N64Co0.5/isoFaces/isoFaces_197.vtk' ]

createImages(isoFacesPath[0],
             isoFacesPath[1],
             os.path.join("S:/Verification/VoFLibary/release/Advection/deformationSphere/Figure" , outputPath))
#S:\Verification\VoFLibary\release\Advection\deformationSphere\plicRDFN\hex\N64Co0.5\isoFaces



#for meshType in os.listdir(algorithm):
#    tmppath = os.path.join(algorithm,meshType)
#    for resolution in os.listdir(tmppath):
#
#        tmppath2 = os.path.join(tmppath,resolution + "/isoFaces")
#        isoFacesPath = list()
#        for isoFaces in os.listdir(tmppath2):
#            finalPath = os.path.join(tmppath2,isoFaces)
#            isoFacesPath.append(finalPath)
#            print finalPath
#        outputPath = meshType + resolution + ".png"
#        print outputPath
#        print os.path.join("S:/Verification/VoFLibary/release/Advection/deformationSphere/Figure" , outputPath)
#        createImages(isoFacesPath[0],
#                     isoFacesPath[1],
#                     os.path.join("S:/Verification/VoFLibary/release/Advection/deformationSphere/Figure" , outputPath))
#                    os.path.join("D:/Paper/VoFLibary/Figure" , outputPath))
    
    
#createImages('./plicRDFN/hex/N128Co0.5/isoFaces/isoFaces_844.vtk',
#             './plicRDFN/hex/N128Co0.5/isoFaces/isoFaces_1683.vtk',
#             './plicRDFN/hex/N128Co0.5/isoFaces/isoFaces_1683.vtk')

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
