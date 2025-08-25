# state file generated using paraview version 5.10.0

# uncomment the following three lines to ensure this script works in future versions
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

import glob
import os
import sys

assert len(sys.argv) == 3

sim_dir = sys.argv[1]
scene_dir = sys.argv[2]

resolution = [1920, 1080]


#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = resolution
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [540.4432220458984, 489.90964698791504, 23.292206287384033]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [588.3569351615023, 489.29887877677095, 1811.3753781352136]
renderView1.CameraFocalPoint = [588.3569351615023, 489.29887877677095, -100.27031399900704]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 724.393214696
renderView1.UseColorPaletteForBackground = 0
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
layout1.SetSize(*resolution)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

filenames = sorted(glob.glob(os.path.join(sim_dir, 'ply', 'ABF_*.ply')))
aBF_000 = PLYReader(registrationName='ABF_000*', FileNames=filenames)

filenames = sorted(glob.glob(os.path.join(sim_dir, 'ply', 'rbc_*.ply')))
rbc_000 = PLYReader(registrationName='rbc_000*', FileNames=filenames)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from aBF_000
aBF_000Display = Show(aBF_000, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
aBF_000Display.Representation = 'Surface'
aBF_000Display.AmbientColor = [0.9098039215686274, 1.0, 0.8274509803921568]
aBF_000Display.ColorArrayName = [None, '']
aBF_000Display.DiffuseColor = [0.9098039215686274, 1.0, 0.8274509803921568]
aBF_000Display.Interpolation = 'PBR'
aBF_000Display.Luminosity = 50.0
aBF_000Display.Diffuse = 0.5
aBF_000Display.Roughness = 0.84
aBF_000Display.SelectTCoordArray = 'None'
aBF_000Display.SelectNormalArray = 'None'
aBF_000Display.SelectTangentArray = 'None'
aBF_000Display.OSPRayScaleFunction = 'PiecewiseFunction'
aBF_000Display.SelectOrientationVectors = 'None'
aBF_000Display.ScaleFactor = 2.2003692626953124
aBF_000Display.SelectScaleArray = 'None'
aBF_000Display.GlyphType = 'Arrow'
aBF_000Display.GlyphTableIndexArray = 'None'
aBF_000Display.GaussianRadius = 0.11001846313476563
aBF_000Display.SetScaleArray = [None, '']
aBF_000Display.ScaleTransferFunction = 'PiecewiseFunction'
aBF_000Display.OpacityArray = [None, '']
aBF_000Display.OpacityTransferFunction = 'PiecewiseFunction'
aBF_000Display.DataAxesGrid = 'GridAxesRepresentation'
aBF_000Display.PolarAxes = 'PolarAxesRepresentation'

# show data from rbc_000
rbc_000Display = Show(rbc_000, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
rbc_000Display.Representation = 'Surface'
rbc_000Display.AmbientColor = [1.0, 0.0, 0.0]
rbc_000Display.ColorArrayName = [None, '']
rbc_000Display.DiffuseColor = [1.0, 0.0, 0.0]
rbc_000Display.Interpolation = 'PBR'
rbc_000Display.Diffuse = 0.5
rbc_000Display.Roughness = 0.84
rbc_000Display.SelectTCoordArray = 'None'
rbc_000Display.SelectNormalArray = 'None'
rbc_000Display.SelectTangentArray = 'None'
rbc_000Display.OSPRayScaleFunction = 'PiecewiseFunction'
rbc_000Display.SelectOrientationVectors = 'None'
rbc_000Display.ScaleFactor = 52.950415039062506
rbc_000Display.SelectScaleArray = 'None'
rbc_000Display.GlyphType = 'Arrow'
rbc_000Display.GlyphTableIndexArray = 'None'
rbc_000Display.GaussianRadius = 2.647520751953125
rbc_000Display.SetScaleArray = [None, '']
rbc_000Display.ScaleTransferFunction = 'PiecewiseFunction'
rbc_000Display.OpacityArray = [None, '']
rbc_000Display.OpacityTransferFunction = 'PiecewiseFunction'
rbc_000Display.DataAxesGrid = 'GridAxesRepresentation'
rbc_000Display.PolarAxes = 'PolarAxesRepresentation'

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
pNG1.Trigger = 'TimeStep'

# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'RenderView1_{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = resolution
pNG1.Writer.Format = 'PNG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(pNG1)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')
