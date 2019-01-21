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
    
    #renderView1.ResetCamera()
    isoFaces_1vtkvtkDisplay = Show(isoFaces_1vtk, renderView1)
    isoFaces_2vtkvtkDisplay = Show(isoFaces_2vtk, renderView1)
    

   
    # show data in view
   
    renderView1.ResetCamera()
    #change interaction mode for render view
    
    # current camera placement for renderView1
    # current camera placement for renderView1
    renderView1.CameraPosition = [0.40737003460526466, 0.5281376093626022, 2.8328371376901633]
    renderView1.CameraFocalPoint = [0.40737003460526466, 0.5281376093626022, 0.45497142150998116]
    renderView1.CameraParallelScale = 0.6154369340437763
    
    SaveScreenshot(outputPath, magnification=1, quality=100, view=renderView1)
        
    # destroy legacyVTKReader2
    Delete(isoFaces_1vtk)
    del isoFaces_1vtk
    
    # destroy legacyVTKReader1
    Delete(isoFaces_2vtk)
    del isoFaces_2vtk
    


os.chdir('S:/Verification/VoFLibary/release/Advection/deformationSphere/')
print os.listdir('./plicRDFN')
algorithm =  './plicRDFN'

#isoFacesPath = [ 'S:/Verification\VoFLibary/release/Advection/deformationSphere/plicRDFN/hex/N64Co0.5/isoFaces/isoFaces_102.vtk','S:/Verification/VoFLibary/release/Advection/deformationSphere/plicRDFN/hex/N64Co0.5/isoFaces/isoFaces_197.vtk' ]
#
#createImages(isoFacesPath[0],
#             isoFacesPath[1],
#             os.path.join("S:/Verification/VoFLibary/release/Advection/deformationSphere/Figure" , outputPath))
#S:\Verification\VoFLibary\release\Advection\deformationSphere\plicRDFN\hex\N64Co0.5\isoFaces



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
        print os.path.join("S:/Verification/VoFLibary/release/Advection/deformationSphere/Figure" , outputPath)
        createImages(isoFacesPath[0],
                     isoFacesPath[1],
                     os.path.join("S:/Verification/VoFLibary/release/Advection/deformationSphere/Figure" , outputPath))


