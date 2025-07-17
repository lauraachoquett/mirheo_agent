# state file generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

import glob
import os
import sys

if len(sys.argv) == 1: # dev mode
    basedir = '/Users/amoudruz/mounts/barry/mirRun/microswimmers/blood_ABF/tube_single_CTC'
    dst = 'extracts_test'
else: # pvbatch mode
    assert len(sys.argv) == 3
    basedir = sys.argv[1]
    dst =  sys.argv[2]

nframes=len(glob.glob(os.path.join(basedir, 'ply', 'rbc*.ply')))

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920, 1080]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [75.0, 95.0, 95.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [75.1761137456564, 15.807764834535789, 107.02740027424898]
renderView1.CameraFocalPoint = [75.1761137456564, 15.807764834535789, 14.899467468261719]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 194.10048943781672
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.EnableRayTracing = 1
renderView1.BackEnd = 'OSPRay pathtracer'
renderView1.SamplesPerPixel = 4
renderView1.EnvironmentalBG = [1.0, 1.0, 1.0]
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1920, 1080)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PLY Reader'
aBF_00 = PLYReader(registrationName='ABF_00*', FileNames=[os.path.join(basedir, 'ply', f'ABF_{i:05d}.ply') for i in range(nframes)])

# create a new 'Plane'
plane2 = Plane(registrationName='Plane2')
plane2.Origin = [-500.0, 0.0, -500.0]
plane2.Point1 = [500.0, 0.0, -500.0]
plane2.Point2 = [-500.0, 0.0, 500.0]

# create a new 'Generate Surface Normals'
generateSurfaceNormals2 = GenerateSurfaceNormals(registrationName='GenerateSurfaceNormals2', Input=aBF_00)

# create a new 'Plane'
plane1 = Plane(registrationName='Plane1')
plane1.Origin = [-500.0, -500.0, 0.0]
plane1.Point1 = [500.0, -500.0, 0.0]
plane1.Point2 = [-500.0, 500.0, 0.0]

# create a new 'PLY Reader'
rbc_00 = PLYReader(registrationName='rbc_00*', FileNames=[os.path.join(basedir, 'ply', f'rbc_{i:05d}.ply') for i in range(nframes)])

# create a new 'Generate Surface Normals'
generateSurfaceNormals1 = GenerateSurfaceNormals(registrationName='GenerateSurfaceNormals1', Input=rbc_00)

# create a new 'Transform'
transform1 = Transform(registrationName='Transform1', Input=generateSurfaceNormals1)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Translate = [149.700592, 0.0, 0.0]


cTC_000 = PLYReader(registrationName='CTC_000*', FileNames=[os.path.join(basedir, 'ply', f'CTC_{i:05d}.ply') for i in range(nframes)])

# create a new 'Generate Surface Normals'
generateSurfaceNormals3 = GenerateSurfaceNormals(registrationName='GenerateSurfaceNormals3', Input=cTC_000)


# create a new 'Cylinder'
cylinder1 = Cylinder(registrationName='Cylinder1')
cylinder1.Resolution = 256
cylinder1.Height = 500.0
cylinder1.Radius = 15.1
cylinder1.Capping = 0

# create a new 'Transform'
transform2 = Transform(registrationName='Transform2', Input=cylinder1)
transform2.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform2.Transform.Translate = [150.0, 16.9700605, 16.9700605]
transform2.Transform.Rotate = [0.0, 0.0, 90.0]

# create a new 'Transform'
transform3 = Transform(registrationName='Transform3', Input=cylinder1)
transform3.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform3.Transform.Translate = [150.0, 16.9700605, 16.9700605]
transform3.Transform.Rotate = [0.0, 0.0, 90.0]
transform3.Transform.Scale = [1.0, 1.02, 1.02]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from generateSurfaceNormals1
generateSurfaceNormals1Display = Show(generateSurfaceNormals1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
generateSurfaceNormals1Display.Representation = 'Surface'
generateSurfaceNormals1Display.AmbientColor = [1.0, 0.0, 0.0]
generateSurfaceNormals1Display.ColorArrayName = [None, '']
generateSurfaceNormals1Display.DiffuseColor = [1.0, 0.0, 0.0]
generateSurfaceNormals1Display.Interpolation = 'PBR'
generateSurfaceNormals1Display.Diffuse = 0.5
generateSurfaceNormals1Display.Roughness = 0.84
generateSurfaceNormals1Display.SelectTCoordArray = 'None'
generateSurfaceNormals1Display.SelectNormalArray = 'Normals'
generateSurfaceNormals1Display.SelectTangentArray = 'None'
generateSurfaceNormals1Display.OSPRayScaleArray = 'Normals'
generateSurfaceNormals1Display.OSPRayScaleFunction = 'PiecewiseFunction'
generateSurfaceNormals1Display.SelectOrientationVectors = 'None'
generateSurfaceNormals1Display.ScaleFactor = 15.21916732788086
generateSurfaceNormals1Display.SelectScaleArray = 'None'
generateSurfaceNormals1Display.GlyphType = 'Arrow'
generateSurfaceNormals1Display.GlyphTableIndexArray = 'None'
generateSurfaceNormals1Display.GaussianRadius = 0.760958366394043
generateSurfaceNormals1Display.SetScaleArray = ['POINTS', 'Normals']
generateSurfaceNormals1Display.ScaleTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals1Display.OpacityArray = ['POINTS', 'Normals']
generateSurfaceNormals1Display.OpacityTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals1Display.DataAxesGrid = 'GridAxesRepresentation'
generateSurfaceNormals1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
generateSurfaceNormals1Display.ScaleTransferFunction.Points = [-0.9999799728393555, 0.0, 0.5, 0.0, 0.999993622303009, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
generateSurfaceNormals1Display.OpacityTransferFunction.Points = [-0.9999799728393555, 0.0, 0.5, 0.0, 0.999993622303009, 1.0, 0.5, 0.0]

# show data from transform1
transform1Display = Show(transform1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
transform1Display.Representation = 'Surface'
transform1Display.AmbientColor = [1.0, 0.0, 0.0]
transform1Display.ColorArrayName = [None, '']
transform1Display.DiffuseColor = [1.0, 0.0, 0.0]
transform1Display.Interpolation = 'PBR'
transform1Display.Diffuse = 0.5
transform1Display.Roughness = 0.84
transform1Display.SelectTCoordArray = 'None'
transform1Display.SelectNormalArray = 'Normals'
transform1Display.SelectTangentArray = 'None'
transform1Display.OSPRayScaleArray = 'Normals'
transform1Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform1Display.SelectOrientationVectors = 'None'
transform1Display.ScaleFactor = 15.219168090820313
transform1Display.SelectScaleArray = 'None'
transform1Display.GlyphType = 'Arrow'
transform1Display.GlyphTableIndexArray = 'None'
transform1Display.GaussianRadius = 0.7609584045410156
transform1Display.SetScaleArray = ['POINTS', 'Normals']
transform1Display.ScaleTransferFunction = 'PiecewiseFunction'
transform1Display.OpacityArray = ['POINTS', 'Normals']
transform1Display.OpacityTransferFunction = 'PiecewiseFunction'
transform1Display.DataAxesGrid = 'GridAxesRepresentation'
transform1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform1Display.ScaleTransferFunction.Points = [-0.9999799728393555, 0.0, 0.5, 0.0, 0.999993622303009, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform1Display.OpacityTransferFunction.Points = [-0.9999799728393555, 0.0, 0.5, 0.0, 0.999993622303009, 1.0, 0.5, 0.0]

# show data from generateSurfaceNormals2
generateSurfaceNormals2Display = Show(generateSurfaceNormals2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
generateSurfaceNormals2Display.Representation = 'Surface'
generateSurfaceNormals2Display.AmbientColor = [1.0, 0.9764705882352941, 0.8117647058823529]
generateSurfaceNormals2Display.ColorArrayName = [None, '']
generateSurfaceNormals2Display.DiffuseColor = [1.0, 0.9764705882352941, 0.8117647058823529]
generateSurfaceNormals2Display.Interpolation = 'PBR'
generateSurfaceNormals2Display.Diffuse = 0.7
generateSurfaceNormals2Display.Roughness = 0.9
generateSurfaceNormals2Display.SelectTCoordArray = 'None'
generateSurfaceNormals2Display.SelectNormalArray = 'Normals'
generateSurfaceNormals2Display.SelectTangentArray = 'None'
generateSurfaceNormals2Display.OSPRayScaleArray = 'Normals'
generateSurfaceNormals2Display.OSPRayScaleFunction = 'PiecewiseFunction'
generateSurfaceNormals2Display.SelectOrientationVectors = 'None'
generateSurfaceNormals2Display.ScaleFactor = 2.24550895690918
generateSurfaceNormals2Display.SelectScaleArray = 'None'
generateSurfaceNormals2Display.GlyphType = 'Arrow'
generateSurfaceNormals2Display.GlyphTableIndexArray = 'None'
generateSurfaceNormals2Display.GaussianRadius = 0.11227544784545898
generateSurfaceNormals2Display.SetScaleArray = ['POINTS', 'Normals']
generateSurfaceNormals2Display.ScaleTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals2Display.OpacityArray = ['POINTS', 'Normals']
generateSurfaceNormals2Display.OpacityTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals2Display.DataAxesGrid = 'GridAxesRepresentation'
generateSurfaceNormals2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
generateSurfaceNormals2Display.ScaleTransferFunction.Points = [-0.9976751804351807, 0.0, 0.5, 0.0, 0.9985284805297852, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
generateSurfaceNormals2Display.OpacityTransferFunction.Points = [-0.9976751804351807, 0.0, 0.5, 0.0, 0.9985284805297852, 1.0, 0.5, 0.0]

# show data from transform2
transform2Display = Show(transform2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
transform2Display.Representation = 'Surface'
transform2Display.ColorArrayName = [None, '']
transform2Display.Opacity = 0.18
transform2Display.SelectTCoordArray = 'TCoords'
transform2Display.SelectNormalArray = 'Normals'
transform2Display.SelectTangentArray = 'None'
transform2Display.OSPRayScaleArray = 'Normals'
transform2Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform2Display.OSPRayMaterial = 'Glass_Thick'
transform2Display.SelectOrientationVectors = 'None'
transform2Display.ScaleFactor = 20.0
transform2Display.SelectScaleArray = 'None'
transform2Display.GlyphType = 'Arrow'
transform2Display.GlyphTableIndexArray = 'None'
transform2Display.GaussianRadius = 1.0
transform2Display.SetScaleArray = ['POINTS', 'Normals']
transform2Display.ScaleTransferFunction = 'PiecewiseFunction'
transform2Display.OpacityArray = ['POINTS', 'Normals']
transform2Display.OpacityTransferFunction = 'PiecewiseFunction'
transform2Display.DataAxesGrid = 'GridAxesRepresentation'
transform2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform2Display.ScaleTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform2Display.OpacityTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# show data from transform3
transform3Display = Show(transform3, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
transform3Display.Representation = 'Surface'
transform3Display.ColorArrayName = [None, '']
transform3Display.Opacity = 0.23
transform3Display.SelectTCoordArray = 'TCoords'
transform3Display.SelectNormalArray = 'Normals'
transform3Display.SelectTangentArray = 'None'
transform3Display.OSPRayScaleArray = 'Normals'
transform3Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform3Display.OSPRayMaterial = 'Glass_Thick'
transform3Display.SelectOrientationVectors = 'None'
transform3Display.ScaleFactor = 51.0
transform3Display.SelectScaleArray = 'None'
transform3Display.GlyphType = 'Arrow'
transform3Display.GlyphTableIndexArray = 'None'
transform3Display.GaussianRadius = 2.5500000000000003
transform3Display.SetScaleArray = ['POINTS', 'Normals']
transform3Display.ScaleTransferFunction = 'PiecewiseFunction'
transform3Display.OpacityArray = ['POINTS', 'Normals']
transform3Display.OpacityTransferFunction = 'PiecewiseFunction'
transform3Display.DataAxesGrid = 'GridAxesRepresentation'
transform3Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform3Display.ScaleTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform3Display.OpacityTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# show data from plane1
plane1Display = Show(plane1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
plane1Display.Representation = 'Surface'
plane1Display.AmbientColor = [0.49411764705882355, 0.6274509803921569, 0.803921568627451]
plane1Display.ColorArrayName = [None, '']
plane1Display.DiffuseColor = [0.49411764705882355, 0.6274509803921569, 0.803921568627451]
plane1Display.Interpolation = 'PBR'
plane1Display.Diffuse = 0.7
plane1Display.Roughness = 0.9
plane1Display.SelectTCoordArray = 'TextureCoordinates'
plane1Display.SelectNormalArray = 'Normals'
plane1Display.SelectTangentArray = 'None'
plane1Display.OSPRayScaleArray = 'Normals'
plane1Display.OSPRayScaleFunction = 'PiecewiseFunction'
plane1Display.SelectOrientationVectors = 'None'
plane1Display.ScaleFactor = 100.0
plane1Display.SelectScaleArray = 'None'
plane1Display.GlyphType = 'Arrow'
plane1Display.GlyphTableIndexArray = 'None'
plane1Display.GaussianRadius = 5.0
plane1Display.SetScaleArray = ['POINTS', 'Normals']
plane1Display.ScaleTransferFunction = 'PiecewiseFunction'
plane1Display.OpacityArray = ['POINTS', 'Normals']
plane1Display.OpacityTransferFunction = 'PiecewiseFunction'
plane1Display.DataAxesGrid = 'GridAxesRepresentation'
plane1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
plane1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
plane1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# show data from plane2
plane2Display = Show(plane2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
plane2Display.Representation = 'Surface'
plane2Display.AmbientColor = [0.788235294117647, 0.788235294117647, 0.788235294117647]
plane2Display.ColorArrayName = [None, '']
plane2Display.DiffuseColor = [0.788235294117647, 0.788235294117647, 0.788235294117647]
plane2Display.Interpolation = 'PBR'
plane2Display.Diffuse = 0.7
plane2Display.Roughness = 0.9
plane2Display.SelectTCoordArray = 'TextureCoordinates'
plane2Display.SelectNormalArray = 'Normals'
plane2Display.SelectTangentArray = 'None'
plane2Display.OSPRayScaleArray = 'Normals'
plane2Display.OSPRayScaleFunction = 'PiecewiseFunction'
plane2Display.SelectOrientationVectors = 'None'
plane2Display.ScaleFactor = 100.0
plane2Display.SelectScaleArray = 'None'
plane2Display.GlyphType = 'Arrow'
plane2Display.GlyphTableIndexArray = 'None'
plane2Display.GaussianRadius = 5.0
plane2Display.SetScaleArray = ['POINTS', 'Normals']
plane2Display.ScaleTransferFunction = 'PiecewiseFunction'
plane2Display.OpacityArray = ['POINTS', 'Normals']
plane2Display.OpacityTransferFunction = 'PiecewiseFunction'
plane2Display.DataAxesGrid = 'GridAxesRepresentation'
plane2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
plane2Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
plane2Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]


# show data from generateSurfaceNormals3
generateSurfaceNormals3Display = Show(generateSurfaceNormals3, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
generateSurfaceNormals3Display.Representation = 'Surface'
generateSurfaceNormals3Display.AmbientColor = [0.6313725490196078, 0.7725490196078432, 0.4745098039215686]
generateSurfaceNormals3Display.ColorArrayName = [None, '']
generateSurfaceNormals3Display.DiffuseColor = [0.6313725490196078, 0.7725490196078432, 0.4745098039215686]
generateSurfaceNormals3Display.Interpolation = 'PBR'
generateSurfaceNormals3Display.Diffuse = 0.5
generateSurfaceNormals3Display.Roughness = 0.7
generateSurfaceNormals3Display.SelectTCoordArray = 'None'
generateSurfaceNormals3Display.SelectNormalArray = 'Normals'
generateSurfaceNormals3Display.SelectTangentArray = 'None'
generateSurfaceNormals3Display.OSPRayScaleArray = 'Normals'
generateSurfaceNormals3Display.OSPRayScaleFunction = 'PiecewiseFunction'
generateSurfaceNormals3Display.SelectOrientationVectors = 'None'
generateSurfaceNormals3Display.ScaleFactor = 1.3473052978515625
generateSurfaceNormals3Display.SelectScaleArray = 'None'
generateSurfaceNormals3Display.GlyphType = 'Arrow'
generateSurfaceNormals3Display.GlyphTableIndexArray = 'None'
generateSurfaceNormals3Display.GaussianRadius = 0.06736526489257813
generateSurfaceNormals3Display.SetScaleArray = ['POINTS', 'Normals']
generateSurfaceNormals3Display.ScaleTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals3Display.OpacityArray = ['POINTS', 'Normals']
generateSurfaceNormals3Display.OpacityTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals3Display.DataAxesGrid = 'GridAxesRepresentation'
generateSurfaceNormals3Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
generateSurfaceNormals3Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
generateSurfaceNormals3Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'RenderView1_%.6ts%cm.png'
pNG1.Writer.ImageResolution = [1920, 1080]
pNG1.Writer.Format = 'PNG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(pNG1)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory=dst)
