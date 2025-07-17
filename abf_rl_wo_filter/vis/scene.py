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

if len(sys.argv) == 3:
    basedir = sys.argv[1]
    dst =  sys.argv[2]
    nframes=len(glob.glob(os.path.join(basedir, 'ply', 'rbc*.ply')))
else:
    basedir = '/scratch/snx3000/amlucas/blood_ABF/capillary/Ht_0.20_seed_12345/'
    dst = 'extracts'
    nframes=499

# domain extents
Lx, Ly, Lz = (149.70059880239523, 33.94011976047905, 33.94011976047905)
rc = 1
R = (Ly - 4*rc) / 2

Lp = 500 # extents of plane and cylinder

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920, 1080]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [75.0, 95.0, 95.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [105.48409807464922, 19.111446915367665, 126.37426616350632]
renderView1.CameraFocalPoint = [105.48409807464922, 19.111446915367665, 14.899467468261719]
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
ABF_PLY = PLYReader(registrationName='ABF_00*', FileNames=[os.path.join(basedir, 'ply', f'abf_{i:05d}.ply') for i in range(nframes)])

# create a new 'Plane'
plane2 = Plane(registrationName='Plane2')
plane2.Origin = [-Lp, 0.0, -Lp]
plane2.Point1 = [+Lp, 0.0, -Lp]
plane2.Point2 = [-Lp, 0.0, +Lp]

# create a new 'Generate Surface Normals'
ABF_SN = GenerateSurfaceNormals(registrationName='GenerateSurfaceNormals2', Input=ABF_PLY)

# create a new 'Plane'
plane1 = Plane(registrationName='Plane1')
plane1.Origin = [-Lp, -Lp, 0.0]
plane1.Point1 = [+Lp, -Lp, 0.0]
plane1.Point2 = [-Lp, +Lp, 0.0]

# create a new 'PLY Reader'
RBC_PLY = PLYReader(registrationName='rbc_00*', FileNames=[os.path.join(basedir, 'ply', f'rbc_{i:05d}.ply') for i in range(nframes)])

# create a new 'Generate Surface Normals'
RBC_SN = GenerateSurfaceNormals(registrationName='GenerateSurfaceNormals1', Input=RBC_PLY)

# create a new 'Transform'
RBC_shifted = Transform(registrationName='RBC_Shifted', Input=RBC_SN)
RBC_shifted.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
RBC_shifted.Transform.Translate = [Lx, 0.0, 0.0]

# create a new 'Cylinder'
cylinder_base = Cylinder(registrationName='Cylinder_Base')
cylinder_base.Resolution = 256
cylinder_base.Height = Lp
cylinder_base.Radius = R + 0.1 # 0.1 for safety: noRBC should go through
cylinder_base.Capping = 0

# create a new 'Transform'
cylinder_inner = Transform(registrationName='Cylinder_Inner', Input=cylinder_base)
cylinder_inner.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
cylinder_inner.Transform.Translate = [Lx, Ly/2, Lz/2]
cylinder_inner.Transform.Rotate = [0.0, 0.0, 90.0]

# create a new 'Transform'
cylinder_outer = Transform(registrationName='Cylinder_Outer', Input=cylinder_base)
cylinder_outer.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
cylinder_outer.Transform.Translate = [Lx, Ly/2, Lz/2]
cylinder_outer.Transform.Rotate = [0.0, 0.0, 90.0]
cylinder_outer.Transform.Scale = [1.0, 1.02, 1.02]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

for rbc_data in [RBC_SN, RBC_shifted]:
    RBC_display = Show(rbc_data, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    RBC_display.Representation = 'Surface'
    RBC_display.AmbientColor = [1.0, 0.0, 0.0]
    RBC_display.ColorArrayName = [None, '']
    RBC_display.DiffuseColor = [1.0, 0.0, 0.0]
    RBC_display.Interpolation = 'PBR'
    RBC_display.Diffuse = 0.5
    RBC_display.Roughness = 0.84
    RBC_display.SelectTCoordArray = 'None'
    RBC_display.SelectNormalArray = 'Normals'
    RBC_display.SelectTangentArray = 'None'
    RBC_display.OSPRayScaleArray = 'Normals'
    RBC_display.OSPRayScaleFunction = 'PiecewiseFunction'
    RBC_display.SelectOrientationVectors = 'None'
    RBC_display.ScaleFactor = 15.21916732788086
    RBC_display.SelectScaleArray = 'None'
    RBC_display.GlyphType = 'Arrow'
    RBC_display.GlyphTableIndexArray = 'None'
    RBC_display.GaussianRadius = 0.760958366394043
    RBC_display.SetScaleArray = ['POINTS', 'Normals']
    RBC_display.ScaleTransferFunction = 'PiecewiseFunction'
    RBC_display.OpacityArray = ['POINTS', 'Normals']
    RBC_display.OpacityTransferFunction = 'PiecewiseFunction'
    RBC_display.DataAxesGrid = 'GridAxesRepresentation'
    RBC_display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    RBC_display.ScaleTransferFunction.Points = [-0.9999799728393555, 0.0, 0.5, 0.0, 0.999993622303009, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    RBC_display.OpacityTransferFunction.Points = [-0.9999799728393555, 0.0, 0.5, 0.0, 0.999993622303009, 1.0, 0.5, 0.0]


# show data from ABF_SN
ABF_display = Show(ABF_SN, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
ABF_display.Representation = 'Surface'
ABF_display.AmbientColor = [1.0, 0.9764705882352941, 0.8117647058823529]
ABF_display.ColorArrayName = [None, '']
ABF_display.DiffuseColor = [1.0, 0.9764705882352941, 0.8117647058823529]
ABF_display.Interpolation = 'PBR'
ABF_display.Diffuse = 0.7
ABF_display.Roughness = 0.9
ABF_display.SelectTCoordArray = 'None'
ABF_display.SelectNormalArray = 'Normals'
ABF_display.SelectTangentArray = 'None'
ABF_display.OSPRayScaleArray = 'Normals'
ABF_display.OSPRayScaleFunction = 'PiecewiseFunction'
ABF_display.SelectOrientationVectors = 'None'
ABF_display.ScaleFactor = 2.24550895690918
ABF_display.SelectScaleArray = 'None'
ABF_display.GlyphType = 'Arrow'
ABF_display.GlyphTableIndexArray = 'None'
ABF_display.GaussianRadius = 0.11227544784545898
ABF_display.SetScaleArray = ['POINTS', 'Normals']
ABF_display.ScaleTransferFunction = 'PiecewiseFunction'
ABF_display.OpacityArray = ['POINTS', 'Normals']
ABF_display.OpacityTransferFunction = 'PiecewiseFunction'
ABF_display.DataAxesGrid = 'GridAxesRepresentation'
ABF_display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
ABF_display.ScaleTransferFunction.Points = [-0.9976751804351807, 0.0, 0.5, 0.0, 0.9985284805297852, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
ABF_display.OpacityTransferFunction.Points = [-0.9976751804351807, 0.0, 0.5, 0.0, 0.9985284805297852, 1.0, 0.5, 0.0]

# show data from cylinder_inner
for cylinder_data in [cylinder_inner, cylinder_outer]:
    cylinder_display = Show(cylinder_data, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    cylinder_display.Representation = 'Surface'
    cylinder_display.ColorArrayName = [None, '']
    cylinder_display.Opacity = 0.18
    cylinder_display.SelectTCoordArray = 'TCoords'
    cylinder_display.SelectNormalArray = 'Normals'
    cylinder_display.SelectTangentArray = 'None'
    cylinder_display.OSPRayScaleArray = 'Normals'
    cylinder_display.OSPRayScaleFunction = 'PiecewiseFunction'
    cylinder_display.OSPRayMaterial = 'Glass_Thick'
    cylinder_display.SelectOrientationVectors = 'None'
    cylinder_display.ScaleFactor = 20.0
    cylinder_display.SelectScaleArray = 'None'
    cylinder_display.GlyphType = 'Arrow'
    cylinder_display.GlyphTableIndexArray = 'None'
    cylinder_display.GaussianRadius = 1.0
    cylinder_display.SetScaleArray = ['POINTS', 'Normals']
    cylinder_display.ScaleTransferFunction = 'PiecewiseFunction'
    cylinder_display.OpacityArray = ['POINTS', 'Normals']
    cylinder_display.OpacityTransferFunction = 'PiecewiseFunction'
    cylinder_display.DataAxesGrid = 'GridAxesRepresentation'
    cylinder_display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    cylinder_display.ScaleTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    cylinder_display.OpacityTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]


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
