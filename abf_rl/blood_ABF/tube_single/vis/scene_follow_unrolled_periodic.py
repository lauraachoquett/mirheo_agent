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
import pandas as pd
import numpy as np
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

# read the ABF coordinates to find out when it was periodically jumping
df = pd.read_csv(os.path.join(basedir, 'obj_stats', 'abf.csv'))
ABF_x = df['comx'].to_numpy()

dx = np.diff(ABF_x)
idx_pos = np.argwhere(dx > +Lx / 2).flatten()
idx_neg = np.argwhere(dx < -Lx / 2).flatten()

periodic_shift = np.zeros_like(ABF_x)
periodic_shift[idx_neg] += Lx
periodic_shift[idx_pos] -= Lx
periodic_shift = np.cumsum(periodic_shift)
print(periodic_shift.tolist())

Lp = 2000 # extents of plane and cylinder

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

# read RBC mesh
ABF_PLY = PLYReader(registrationName='ABF_00*', FileNames=[os.path.join(basedir, 'ply', f'abf_{i:05d}.ply') for i in range(nframes)])
ABF_SN = GenerateSurfaceNormals(registrationName='GenerateSurfaceNormals2', Input=ABF_PLY)

# wrap the ABF into a transform so we can later shift its position according to the periodicity
ABF_unrolled = Transform(registrationName='ABF_unrolled', Input=ABF_SN)
ABF_unrolled.Transform = 'Transform'
ABF_unrolled.Transform.Translate = [0.0, 0.0, 0.0]

# read RBC mesh
RBC_PLY = PLYReader(registrationName='rbc_00*', FileNames=[os.path.join(basedir, 'ply', f'rbc_{i:05d}.ply') for i in range(nframes)])
RBC_SN = GenerateSurfaceNormals(registrationName='GenerateSurfaceNormals1', Input=RBC_PLY)

RBC_unrolled = Transform(registrationName='RBC_unrolled', Input=RBC_SN)
RBC_unrolled.Transform = 'Transform'
RBC_unrolled.Transform.Translate = [0.0, 0.0, 0.0]

RBC_unrolled_next = Transform(registrationName='RBC_unrolled_next', Input=RBC_unrolled)
RBC_unrolled_next.Transform = 'Transform'
RBC_unrolled_next.Transform.Translate = [Lx, 0.0, 0.0]

RBC_unrolled_prev = Transform(registrationName='RBC_unrolled_prev', Input=RBC_unrolled)
RBC_unrolled_prev.Transform = 'Transform'
RBC_unrolled_prev.Transform.Translate = [-Lx, 0.0, 0.0]


# create a new 'Plane'
plane_back = Plane(registrationName='Plane_Back')
plane_back.Origin = [-Lp, -Lp, 0.0]
plane_back.Point1 = [+Lp, -Lp, 0.0]
plane_back.Point2 = [-Lp, +Lp, 0.0]

# create a new 'Plane'
plane_bottom = Plane(registrationName='Plane_Bottom')
plane_bottom.Origin = [-Lp, 0.0, -Lp]
plane_bottom.Point1 = [+Lp, 0.0, -Lp]
plane_bottom.Point2 = [-Lp, 0.0, +Lp]
plane_bottom.XResolution = 4*2048
plane_bottom.YResolution = 4*2048


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

for rbc_data in [RBC_unrolled, RBC_unrolled_prev, RBC_unrolled_next]:
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


ABF_display = Show(ABF_unrolled, renderView1, 'GeometryRepresentation')

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


# show data from plane_back
plane_back_display = Show(plane_back, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
plane_back_display.Representation = 'Surface'
plane_back_display.AmbientColor = [0.49411764705882355, 0.6274509803921569, 0.803921568627451]
plane_back_display.ColorArrayName = [None, '']
plane_back_display.DiffuseColor = [0.49411764705882355, 0.6274509803921569, 0.803921568627451]
plane_back_display.Interpolation = 'PBR'
plane_back_display.Diffuse = 0.7
plane_back_display.Roughness = 0.9
plane_back_display.SelectTCoordArray = 'TextureCoordinates'
plane_back_display.SelectNormalArray = 'Normals'
plane_back_display.SelectTangentArray = 'None'
plane_back_display.OSPRayScaleArray = 'Normals'
plane_back_display.OSPRayScaleFunction = 'PiecewiseFunction'
plane_back_display.SelectOrientationVectors = 'None'
plane_back_display.ScaleFactor = 100.0
plane_back_display.SelectScaleArray = 'None'
plane_back_display.GlyphType = 'Arrow'
plane_back_display.GlyphTableIndexArray = 'None'
plane_back_display.GaussianRadius = 5.0
plane_back_display.SetScaleArray = ['POINTS', 'Normals']
plane_back_display.ScaleTransferFunction = 'PiecewiseFunction'
plane_back_display.OpacityArray = ['POINTS', 'Normals']
plane_back_display.OpacityTransferFunction = 'PiecewiseFunction'
plane_back_display.DataAxesGrid = 'GridAxesRepresentation'
plane_back_display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
plane_back_display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
plane_back_display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# create chess pattern on the bottom plane
chess_L = Lx / 4
chess_k = (2*Lp) / chess_L * 2 * np.pi

plane_bottom_with_pattern = Calculator(registrationName='plane_bottom_with_pattern', Input=plane_bottom)
plane_bottom_with_pattern.ResultTCoords = 1
plane_bottom_with_pattern.ResultArrayName = 'pattern'
plane_bottom_with_pattern.Function = f'sin({chess_k} * TextureCoordinates_X) * sin({chess_k} * TextureCoordinates_Y)'

plane_bottom_display = Show(plane_bottom_with_pattern, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'pattern'
patternLUT = GetColorTransferFunction('pattern')
patternLUT.RGBPoints = [-1.0, 0.8, 0.8, 0.8,
                        0.0, 0.8, 0.8, 0.8,
                        1.0, 0.7, 0.7, 0.7]
patternLUT.ColorSpace = 'Step'
patternLUT.NanColor = [0.803922, 0.0, 0.803922]
patternLUT.ScalarRangeInitialized = 1.0

plane_bottom_display.Representation = 'Surface'
plane_bottom_display.ColorArrayName = ['POINTS', 'pattern']
plane_bottom_display.LookupTable = patternLUT
plane_bottom_display.Interpolation = 'PBR'
plane_bottom_display.Diffuse = 0.7
plane_bottom_display.Roughness = 0.9
plane_bottom_display.SelectTCoordArray = 'TextureCoordinates'
plane_bottom_display.SelectNormalArray = 'Normals'
plane_bottom_display.SelectTangentArray = 'None'
plane_bottom_display.OSPRayScaleArray = 'Normals'
plane_bottom_display.OSPRayScaleFunction = 'PiecewiseFunction'
plane_bottom_display.SelectOrientationVectors = 'None'
plane_bottom_display.ScaleFactor = 100.0
plane_bottom_display.SelectScaleArray = 'None'
plane_bottom_display.GlyphType = 'Arrow'
plane_bottom_display.GlyphTableIndexArray = 'None'
plane_bottom_display.GaussianRadius = 5.0
plane_bottom_display.SetScaleArray = ['POINTS', 'Normals']
plane_bottom_display.ScaleTransferFunction = 'PiecewiseFunction'
plane_bottom_display.OpacityArray = ['POINTS', 'Normals']
plane_bottom_display.OpacityTransferFunction = 'PiecewiseFunction'
plane_bottom_display.DataAxesGrid = 'GridAxesRepresentation'
plane_bottom_display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
plane_bottom_display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
plane_bottom_display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

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


scene = GetAnimationScene()
scene.NumberOfFrames = nframes
scene.PlayMode =  'Sequence'

# Adding now a PythonAnimationCue()

pyQ01 = PythonAnimationCue()
pyQ01.Script=f"""from paraview.simple import *

positions = {ABF_x.tolist()}
jumps = {periodic_shift.tolist()}

def start_cue(self): pass

def tick(self):
    tk = GetTimeKeeper()
    t = int(tk.Time)
    shiftx = jumps[t]
    abf = FindSource('ABF_unrolled')
    abf.Transform.Translate[0] = shiftx
    rbc = FindSource('RBC_unrolled')
    rbc.Transform.Translate[0] = shiftx
    abf_x = positions[t] + shiftx
    view = GetRenderView()
    view.CameraPosition[0] = abf_x
    view.CameraFocalPoint[0] = abf_x


def end_cue(self): pass"""
pyQ01.StartTime=0.0
pyQ01.EndTime=1.0
scene.Cues.append(pyQ01)


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory=dst)
