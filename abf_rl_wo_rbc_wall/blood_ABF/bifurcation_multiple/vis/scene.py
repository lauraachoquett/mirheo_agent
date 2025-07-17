# state file generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

import glob
import os
import sys

if len(sys.argv) == 3:
    basedir = sys.argv[1]
    out = sys.argv[2]
else:
    basedir = '/scratch/snx3000/amlucas/blood_ABF/swarm_pipe/Re_0.5_Ht_0.20_L_10_P_4_V_5_mm_per_s_IC_rnd,300,123112_seed_123112_novisc'
    out = 'extracts'


res = [1920, 1080]

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = res
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [259.61173248291016, 89.59198951721191, 40.87903308868408]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [259.61173248291016, 89.59198951721191, 456.4499685816073]
renderView1.CameraFocalPoint = [259.61173248291016, 89.59198951721191, 40.87903308868408]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 278.97690258281295
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.EnableRayTracing = 1
renderView1.BackEnd = 'OSPRay pathtracer'
renderView1.SamplesPerPixel = 4
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(*res)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Plane'
plane1 = Plane(registrationName='Plane1')
plane1.Origin = [-5000.0, -5000.0, -200.0]
plane1.Point1 = [5000.0, -5000.0, -200.0]
plane1.Point2 = [-5000.0, 5000.0, -200.0]

rbc_fnames = sorted(glob.glob(os.path.join(basedir, 'ply', 'rbc_*.ply')))
abf_fnames = sorted(glob.glob(os.path.join(basedir, 'ply', 'ABF_*.ply')))

# create a new 'PLY Reader'
rbc_00 = PLYReader(registrationName='rbc_00*', FileNames=rbc_fnames)

# create a new 'Generate Surface Normals'
generateSurfaceNormals1 = GenerateSurfaceNormals(registrationName='GenerateSurfaceNormals1', Input=rbc_00)

# create a new 'PLY Reader'
aBF_00 = PLYReader(registrationName='ABF_00*', FileNames=abf_fnames)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from aBF_00
aBF_00Display = Show(aBF_00, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
aBF_00Display.Representation = 'Surface'
aBF_00Display.AmbientColor = [1.0, 0.9764705882352941, 0.8117647058823529]
aBF_00Display.ColorArrayName = [None, '']
aBF_00Display.DiffuseColor = [1.0, 0.9764705882352941, 0.8117647058823529]
aBF_00Display.Interpolation = 'PBR'
aBF_00Display.Diffuse = 0.7
aBF_00Display.Roughness = 0.9
aBF_00Display.SelectTCoordArray = 'None'
aBF_00Display.SelectNormalArray = 'None'
aBF_00Display.SelectTangentArray = 'None'
aBF_00Display.OSPRayScaleFunction = 'PiecewiseFunction'
aBF_00Display.SelectOrientationVectors = 'None'
aBF_00Display.ScaleFactor = 17.118035888671876
aBF_00Display.SelectScaleArray = 'None'
aBF_00Display.GlyphType = 'Arrow'
aBF_00Display.GlyphTableIndexArray = 'None'
aBF_00Display.GaussianRadius = 0.8559017944335938
aBF_00Display.SetScaleArray = [None, '']
aBF_00Display.ScaleTransferFunction = 'PiecewiseFunction'
aBF_00Display.OpacityArray = [None, '']
aBF_00Display.OpacityTransferFunction = 'PiecewiseFunction'
aBF_00Display.DataAxesGrid = 'GridAxesRepresentation'
aBF_00Display.PolarAxes = 'PolarAxesRepresentation'

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
generateSurfaceNormals1Display.ScaleFactor = 52.54725799560547
generateSurfaceNormals1Display.SelectScaleArray = 'None'
generateSurfaceNormals1Display.GlyphType = 'Arrow'
generateSurfaceNormals1Display.GlyphTableIndexArray = 'None'
generateSurfaceNormals1Display.GaussianRadius = 2.6273628997802736
generateSurfaceNormals1Display.SetScaleArray = ['POINTS', 'Normals']
generateSurfaceNormals1Display.ScaleTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals1Display.OpacityArray = ['POINTS', 'Normals']
generateSurfaceNormals1Display.OpacityTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals1Display.DataAxesGrid = 'GridAxesRepresentation'
generateSurfaceNormals1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
generateSurfaceNormals1Display.ScaleTransferFunction.Points = [-0.9999998211860657, 0.0, 0.5, 0.0, 0.9999992847442627, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
generateSurfaceNormals1Display.OpacityTransferFunction.Points = [-0.9999998211860657, 0.0, 0.5, 0.0, 0.9999992847442627, 1.0, 0.5, 0.0]

# show data from plane1
plane1Display = Show(plane1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
plane1Display.Representation = 'Surface'
plane1Display.AmbientColor = [0.0, 0.0, 0.0]
plane1Display.ColorArrayName = [None, '']
plane1Display.DiffuseColor = [0.0, 0.0, 0.0]
plane1Display.SelectTCoordArray = 'TextureCoordinates'
plane1Display.SelectNormalArray = 'Normals'
plane1Display.SelectTangentArray = 'None'
plane1Display.OSPRayScaleArray = 'Normals'
plane1Display.OSPRayScaleFunction = 'PiecewiseFunction'
plane1Display.SelectOrientationVectors = 'None'
plane1Display.ScaleFactor = 1000.0
plane1Display.SelectScaleArray = 'None'
plane1Display.GlyphType = 'Arrow'
plane1Display.GlyphTableIndexArray = 'None'
plane1Display.GaussianRadius = 50.0
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

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'RenderView1_%.6ts%cm.png'
pNG1.Writer.ImageResolution = res
pNG1.Writer.Format = 'PNG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(pNG1)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory=out)
