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
    
    # create a new 'Slice'
    slice1 = Slice(Input=isoFaces_1vtk)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    
    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = [0.4673001244664192, 0.0, 0.46639133244752884]
    slice1.SliceType.Normal = [0.0, 1.0, 0.0]
    
    # create a new 'Tube'
    tube1 = Tube(Input=slice1)
    tube1.Scalars = [None, '']
    tube1.Vectors = [None, '']
    tube1.Radius = 0.001
    tube1Display = Show(tube1, renderView1)
    tube1Display.DiffuseColor = [0.0, 0.0, 0.0]
    
    #Hide(isoFaces_844vtk, renderView1)
    #
    ## change solid color
    #slice1Display.DiffuseColor = [0.0, 0.0, 0.0]
    
    # change representation type
#    slice1Display.SetRepresentationType('Wireframe')
#    slice1Display.LineWidth = 5.5
#    
#    # change solid color
#    slice1Display.AmbientColor = [0.0, 0.0, 0.0]
    
    #
    #
    # create a new 'Slice'
    slice2 = Slice(Input=isoFaces_2vtk)
    slice2.SliceType = 'Plane'
    slice2.SliceOffsetValues = [0.0]
    
    # init the 'Plane' selected for 'SliceType'
    slice2.SliceType.Origin = [0.4673001244664192, 0.0, 0.46639133244752884]
    slice2.SliceType.Normal = [0.0, 1.0, 0.0]
    
        # create a new 'Tube'
    tube2 = Tube(Input=slice2)
    tube2.Scalars = [None, '']
    tube2.Vectors = [None, '']
    tube2.Radius = 0.001
    tube2Display = Show(tube2, renderView1)
    tube2Display.DiffuseColor = [0.0, 0.0, 0.0]
    
    # show data in view
    #slice2Display = Show(slice2, renderView1)
    
    #Hide(isoFaces_844vtk, renderView1)
    #
    ## change solid color
    #slice1Display.DiffuseColor = [0.0, 0.0, 0.0]
    
    # change representation type
#    slice2Display.SetRepresentationType('Wireframe')
#    slice2Display.LineWidth = 5.5
#    
#    # change solid color
#    slice2Display.AmbientColor = [0.0, 0.0, 0.0]
#    
#    # toggle 3D widget visibility (only when running from the GUI)
#    Hide3DWidgets(proxy=slice1.SliceType)
#    Hide3DWidgets(proxy=slice2.SliceType)
    
    #### saving camera placements for all active views
    
    # current camera placement for renderView1
    #renderView1.CameraPosition = [0.4479870978587231, -1.6659052686488176, 0.46639133244752884]
    #renderView1.CameraFocalPoint = [0.4794869806643196, 1.0512152424756371, 0.46639133244752884]
    #renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    #renderView1.CameraParallelScale = 0.7032897925160194
    
    renderView1.ResetCamera()
    #change interaction mode for render view
    renderView1.InteractionMode = '2D'
    
    #### saving camera placements for all active views
    
    # current camera placement for renderView1
    renderView1.CameraPosition = [0.4441313005842177, -1.9984990730751364, 0.49652717262506485]
    renderView1.CameraFocalPoint = [0.4673001244664192, 0.0, 0.49652717262506485]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraParallelScale = 0.517284379802796
    
    SaveScreenshot(outputPath, magnification=1, quality=100, view=renderView1)
    
    # destroy slice1
    Delete(tube1)
    del tube1
    
    # destroy slice1
    Delete(tube2)
    del tube2
    
    # destroy slice1
    Delete(slice2)
    del slice2
    
    # destroy slice1
    Delete(slice1)
    del slice1
    
    # destroy legacyVTKReader2
    Delete(isoFaces_1vtk)
    del isoFaces_1vtk
    
    # destroy legacyVTKReader1
    Delete(isoFaces_2vtk)
    del isoFaces_2vtk
    

#os.chdir('S:\\Verification\\VoFLibary\\release\\Advection\\vortexShearedDisc')
#testVar = raw_input()
#os.chdir('S:\\Verification\\VoFLibary\\release\\Advection\\vortexShearedDisc')

# create a new 'Legacy VTK Reader'
#isoFaces_1vtk = LegacyVTKReader(FileNames=['S:\\Verification\\VoFLibary\\release\\Advection\\vortexShearedDisc\\plicRDFN\\hex\\N128Co0.5\\isoFaces\\isoFaces_844.vtk'])
#isoFaces_2vtk = LegacyVTKReader(FileNames=['S:\\Verification\\VoFLibary\\release\\Advection\\vortexShearedDisc\\plicRDFN\\hex\\N128Co0.5\\isoFaces\\isoFaces_1683.vtk'])
print os.listdir('./plicRDFN')
algorithm =  './plicRDFN'

for meshType in os.listdir(algorithm):
    tmppath = os.path.join(algorithm,meshType)
    for resolution in os.listdir(tmppath):

        tmppath2 = os.path.join(tmppath,resolution + "/isoFaces")
        isoFacesPath = list()
        for isoFaces in os.listdir(tmppath2):
            finalPath = os.path.join(tmppath2,isoFaces)
            isoFacesPath.append(finalPath)
            print finalPath
        outputPath = meshType + resolution + ".png"
        print outputPath
        createImages(isoFacesPath[0],
                     isoFacesPath[1],
                     os.path.join("S:/Verification/VoFLibary/release/Advection/vortexShearedDisc/Figure" , outputPath))
#                    os.path.join("D:/Paper/VoFLibary/Figure" , outputPath))
    
    
#createImages('./plicRDFN/hex/N128Co0.5/isoFaces/isoFaces_844.vtk',
#             './plicRDFN/hex/N128Co0.5/isoFaces/isoFaces_1683.vtk',
#             './plicRDFN/hex/N128Co0.5/isoFaces/isoFaces_1683.vtk')

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
