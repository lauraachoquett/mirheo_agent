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
renderView1.CenterOfRotation = [148.56503295898438, 39.55927228927612, 39.531758308410645]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [148.56503295898438, 39.55927228927612, 331.4688570607515]
renderView1.CameraFocalPoint = [148.56503295898438, 39.55927228927612, 39.531758308410645]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 161.96717208441865
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
plane1.Origin = [-10000.0, -10000.0, -500.0]
plane1.Point1 = [10000.0, -10000.0, -500.0]
plane1.Point2 = [-10000.0, 10000.0, -500.0]


rbc_fnames = sorted(glob.glob(os.path.join(basedir, 'ply', 'rbc_*.ply')))
abf_fnames = sorted(glob.glob(os.path.join(basedir, 'ply', 'ABF_*.ply')))

# create a new 'PLY Reader'
rbc_00 = PLYReader(registrationName='rbc_00*', FileNames=rbc_fnames)

# create a new 'PLY Reader'
aBF_00 = PLYReader(registrationName='ABF_00*', FileNames=abf_fnames)

# create a new 'Generate Surface Normals'
generateSurfaceNormals2 = GenerateSurfaceNormals(registrationName='GenerateSurfaceNormals2', Input=aBF_00)

# create a new 'Generate Surface Normals'
generateSurfaceNormals1 = GenerateSurfaceNormals(registrationName='GenerateSurfaceNormals1', Input=rbc_00)

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
generateSurfaceNormals1Display.ScaleFactor = 30.771215820312502
generateSurfaceNormals1Display.SelectScaleArray = 'None'
generateSurfaceNormals1Display.GlyphType = 'Arrow'
generateSurfaceNormals1Display.GlyphTableIndexArray = 'None'
generateSurfaceNormals1Display.GaussianRadius = 1.538560791015625
generateSurfaceNormals1Display.SetScaleArray = ['POINTS', 'Normals']
generateSurfaceNormals1Display.ScaleTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals1Display.OpacityArray = ['POINTS', 'Normals']
generateSurfaceNormals1Display.OpacityTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals1Display.DataAxesGrid = 'GridAxesRepresentation'
generateSurfaceNormals1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
generateSurfaceNormals1Display.ScaleTransferFunction.Points = [-0.9999959468841553, 0.0, 0.5, 0.0, 0.9999995827674866, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
generateSurfaceNormals1Display.OpacityTransferFunction.Points = [-0.9999959468841553, 0.0, 0.5, 0.0, 0.9999995827674866, 1.0, 0.5, 0.0]

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
generateSurfaceNormals2Display.ScaleFactor = 31.991500854492188
generateSurfaceNormals2Display.SelectScaleArray = 'None'
generateSurfaceNormals2Display.GlyphType = 'Arrow'
generateSurfaceNormals2Display.GlyphTableIndexArray = 'None'
generateSurfaceNormals2Display.GaussianRadius = 1.5995750427246094
generateSurfaceNormals2Display.SetScaleArray = ['POINTS', 'Normals']
generateSurfaceNormals2Display.ScaleTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals2Display.OpacityArray = ['POINTS', 'Normals']
generateSurfaceNormals2Display.OpacityTransferFunction = 'PiecewiseFunction'
generateSurfaceNormals2Display.DataAxesGrid = 'GridAxesRepresentation'
generateSurfaceNormals2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
generateSurfaceNormals2Display.ScaleTransferFunction.Points = [-0.9987668395042419, 0.0, 0.5, 0.0, 0.9986069202423096, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
generateSurfaceNormals2Display.OpacityTransferFunction.Points = [-0.9987668395042419, 0.0, 0.5, 0.0, 0.9986069202423096, 1.0, 0.5, 0.0]



# show data from plane1
plane1Display = Show(plane1, renderView1, 'GeometryRepresentation')
plane_color = [0.0, 0.0, 0.0]

# trace defaults for the display properties.
plane1Display.Representation = 'Surface'
plane1Display.ColorArrayName = [None, '']
plane1Display.SelectTCoordArray = 'TextureCoordinates'
plane1Display.SelectNormalArray = 'Normals'
plane1Display.SelectTangentArray = 'None'
plane1Display.OSPRayScaleArray = 'Normals'
plane1Display.OSPRayScaleFunction = 'PiecewiseFunction'
plane1Display.SelectOrientationVectors = 'None'
plane1Display.ScaleFactor = 2000.0
plane1Display.SelectScaleArray = 'None'
plane1Display.GlyphType = 'Arrow'
plane1Display.GlyphTableIndexArray = 'None'
plane1Display.GaussianRadius = 100.0
plane1Display.SetScaleArray = ['POINTS', 'Normals']
plane1Display.ScaleTransferFunction = 'PiecewiseFunction'
plane1Display.OpacityArray = ['POINTS', 'Normals']
plane1Display.OpacityTransferFunction = 'PiecewiseFunction'
plane1Display.DataAxesGrid = 'GridAxesRepresentation'
plane1Display.PolarAxes = 'PolarAxesRepresentation'
plane1Display.AmbientColor = plane_color
plane1Display.DiffuseColor = plane_color


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
