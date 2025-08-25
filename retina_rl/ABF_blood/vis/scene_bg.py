from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import glob
import os
import sys

assert len(sys.argv) == 4

sim_dir = sys.argv[1]
scene_dir = sys.argv[2]
out_dir = sys.argv[3]

resolution = [1920, 1080]

target_pos0 = [213.0/576.0, 395.0/528.0, 12.0/24.0]
target_radius0 = 5/576.0
domain = [1071.509572395379, 982.2171080290974, 44.646232183140796]

target_pos = [x * L for x, L in zip(target_pos0, domain)]
target_radius = target_radius0 * domain[0]

materialLibrary1 = GetMaterialLibrary()

renderView1 = CreateView('RenderView')
renderView1.ViewSize = resolution
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [1072.823230852701, 516.7811596541596, 13.49709701538086]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [1072.823230852701, 516.7811596541596, 1912.3926247269505]
renderView1.CameraFocalPoint = [1072.823230852701, 516.7811596541596, -2667.712921179577]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1185.4185438602867
renderView1.UseColorPaletteForBackground = 0
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.EnableRayTracing = 1
renderView1.BackEnd = 'OSPRay pathtracer'
renderView1.SamplesPerPixel = 4
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(*resolution)

SetActiveView(renderView1)

# Background image

background_png = PNGSeriesReader(registrationName='background.png', FileNames=[os.path.join(scene_dir, 'background.png')])

background = Transform(registrationName='Background', Input=background_png)
background.Transform = 'Transform'
background.Transform.Translate = [17.373230852701, -22.718840345840363, 0.0]
background.Transform.Scale = [1.1, 1.0, 0.0]

# RBCs

filenames = sorted(glob.glob(os.path.join(sim_dir, 'ply', 'rbc_*.ply')))
ply_rbc = PLYReader(registrationName='rbc_ply', FileNames=filenames)

clip_rbc = Clip(registrationName='Clip_Rbc', Input=ply_rbc)
clip_rbc.ClipType = 'Plane'
clip_rbc.HyperTreeGridClipper = 'Plane'
clip_rbc.Scalars = ['POINTS', '']
clip_rbc.ClipType.Origin = [990.0, 491.35726737976074, 23.296895027160645]
clip_rbc.HyperTreeGridClipper.Origin = [540.159984588623, 491.35726737976074, 23.296895027160645]

rbc = Transform(registrationName='Rbc', Input=clip_rbc)
rbc.Transform = 'Transform'

# ABF

filenames = sorted(glob.glob(os.path.join(sim_dir, 'ply', 'ABF_*.ply')))
ply_abf = PLYReader(registrationName='ABF_ply', FileNames=filenames)

abf = Transform(registrationName='abf', Input=ply_abf)
abf.Transform = 'Transform'

# Target

sphere_target = Sphere(registrationName='sphere_target')
sphere_target.Center = target_pos
sphere_target.Radius = target_radius
sphere_target.ThetaResolution = 256
sphere_target.PhiResolution = 256

target = Transform(registrationName='target', Input=sphere_target)
target.Transform = 'Transform'


for t in [rbc, abf, target]:
    t.Transform.Translate = [620.0, 300.0, 0.0]
    t.Transform.Rotate = [0.0, 0.0, -20.0]
    t.Transform.Scale = [0.73, 0.73, 0.73]



rbc_display = Show(rbc, renderView1, 'UnstructuredGridRepresentation')
rbc_display.Representation = 'Surface'
rbc_display.AmbientColor = [1.0, 0.0, 0.0]
rbc_display.ColorArrayName = [None, '']
rbc_display.DiffuseColor = [1.0, 0.0, 0.0]
rbc_display.Interpolation = 'PBR'
rbc_display.Diffuse = 0.5
rbc_display.Roughness = 0.84
rbc_display.SelectTCoordArray = 'None'
rbc_display.SelectNormalArray = 'None'
rbc_display.SelectTangentArray = 'None'
rbc_display.OSPRayScaleFunction = 'PiecewiseFunction'
rbc_display.SelectOrientationVectors = 'None'
rbc_display.ScaleFactor = 98.89037857055665
rbc_display.SelectScaleArray = 'None'
rbc_display.GlyphType = 'Arrow'
rbc_display.GlyphTableIndexArray = 'None'
rbc_display.GaussianRadius = 4.944518928527832
rbc_display.SetScaleArray = [None, '']
rbc_display.ScaleTransferFunction = 'PiecewiseFunction'
rbc_display.OpacityArray = [None, '']
rbc_display.OpacityTransferFunction = 'PiecewiseFunction'
rbc_display.DataAxesGrid = 'GridAxesRepresentation'
rbc_display.PolarAxes = 'PolarAxesRepresentation'
rbc_display.ScalarOpacityUnitDistance = 7.478953101683261
rbc_display.OpacityArrayName = [None, '']


background_display = Show(background, renderView1, 'StructuredGridRepresentation')
background_display.Representation = 'Surface'
background_display.ColorArrayName = ['POINTS', 'PNGImage']
background_display.MapScalars = 0
background_display.SelectTCoordArray = 'None'
background_display.SelectNormalArray = 'None'
background_display.SelectTangentArray = 'None'
background_display.OSPRayScaleArray = 'PNGImage'
background_display.OSPRayScaleFunction = 'PiecewiseFunction'
background_display.SelectOrientationVectors = 'None'
background_display.ScaleFactor = 191.9
background_display.SelectScaleArray = 'None'
background_display.GlyphType = 'Arrow'
background_display.GlyphTableIndexArray = 'PNGImage'
background_display.GaussianRadius = 9.595
background_display.SetScaleArray = ['POINTS', 'PNGImage']
background_display.ScaleTransferFunction = 'PiecewiseFunction'
background_display.OpacityArray = ['POINTS', 'PNGImage']
background_display.OpacityTransferFunction = 'PiecewiseFunction'
background_display.DataAxesGrid = 'GridAxesRepresentation'
background_display.PolarAxes = 'PolarAxesRepresentation'
background_display.ScalarOpacityUnitDistance = 17.272776998804524
background_display.ScaleTransferFunction.Points = [6.0, 0.0, 0.5, 0.0, 255.0, 1.0, 0.5, 0.0]
background_display.OpacityTransferFunction.Points = [6.0, 0.0, 0.5, 0.0, 255.0, 1.0, 0.5, 0.0]


abf_display = Show(abf, renderView1, 'GeometryRepresentation')
abf_display.Representation = 'Surface'
abf_display.AmbientColor = [0.9098039215686274, 1.0, 0.8274509803921568]
abf_display.ColorArrayName = [None, '']
abf_display.DiffuseColor = [0.9098039215686274, 1.0, 0.8274509803921568]
abf_display.Interpolation = 'PBR'
abf_display.Luminosity = 50.0
abf_display.Diffuse = 0.5
abf_display.Roughness = 0.84
abf_display.SelectTCoordArray = 'None'
abf_display.SelectNormalArray = 'None'
abf_display.SelectTangentArray = 'None'
abf_display.OSPRayScaleFunction = 'PiecewiseFunction'
abf_display.SelectOrientationVectors = 'None'
abf_display.ScaleFactor = 2.2003692626953124
abf_display.SelectScaleArray = 'None'
abf_display.GlyphType = 'Arrow'
abf_display.GlyphTableIndexArray = 'None'
abf_display.GaussianRadius = 0.11001846313476563
abf_display.SetScaleArray = [None, '']
abf_display.ScaleTransferFunction = 'PiecewiseFunction'
abf_display.OpacityArray = [None, '']
abf_display.OpacityTransferFunction = 'PiecewiseFunction'
abf_display.DataAxesGrid = 'GridAxesRepresentation'
abf_display.PolarAxes = 'PolarAxesRepresentation'


target_display = Show(target, renderView1, 'GeometryRepresentation')
target_display.Representation = 'Surface'
target_display.AmbientColor = [0.9098039215686274, 1.0, 0.8274509803921568]
target_display.ColorArrayName = [None, '']
target_display.DiffuseColor = [0.9098039215686274, 1.0, 0.8274509803921568]
target_display.Interpolation = 'PBR'
target_display.Luminosity = 50.0
target_display.Diffuse = 0.5
target_display.Roughness = 0.84
target_display.SelectTCoordArray = 'None'
target_display.SelectNormalArray = 'None'
target_display.SelectTangentArray = 'None'
target_display.OSPRayScaleFunction = 'PiecewiseFunction'
target_display.SelectOrientationVectors = 'None'
target_display.ScaleFactor = 2.2003692626953124
target_display.SelectScaleArray = 'None'
target_display.GlyphType = 'Arrow'
target_display.GlyphTableIndexArray = 'None'
target_display.GaussianRadius = 0.11001846313476563
target_display.SetScaleArray = [None, '']
target_display.ScaleTransferFunction = 'PiecewiseFunction'
target_display.OpacityArray = [None, '']
target_display.OpacityTransferFunction = 'PiecewiseFunction'
target_display.DataAxesGrid = 'GridAxesRepresentation'
target_display.PolarAxes = 'PolarAxesRepresentation'


pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
pNG1.Trigger = 'TimeStep'
pNG1.Writer.FileName = 'RenderView1_{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = resolution
pNG1.Writer.Format = 'PNG'

SetActiveSource(pNG1)

if __name__ == '__main__':
    SaveExtracts(ExtractsOutputDirectory=out_dir)
