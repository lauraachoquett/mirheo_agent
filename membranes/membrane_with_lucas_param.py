#!/usr/bin/env python

import mirheo as mir
from dpdprops import KantorWLCRBCDefaultParams
import pint
import trimesh
import numpy as np


def rescale_by_area(mesh, A):
    mesh.vertices *= np.sqrt(A / mesh.area)
    return mesh


ureg = pint.UnitRegistry()


params = KantorWLCRBCDefaultParams(
    ureg, mu=4.058 * ureg.uN / ureg.m, kappab=5.906e-19 * ureg.J, x0=0.416
)


dt = 0.001

ranks = (1, 1, 1)
domain = (12, 12, 12)

u = mir.Mirheo(ranks, domain, debug_level=3, log_filename="log")

# we need to first create a mesh before initializing the membrane vector

mesh_ini = trimesh.load_mesh("data/rbc_mesh.off", process=False)
RA = 1
A0 = 4 * np.pi * RA**2
mesh_ini = rescale_by_area(mesh_ini, A0)

mesh_rbc = mir.ParticleVectors.MembraneMesh(
    mesh_ini.vertices.tolist(), mesh_ini.faces.tolist()
)


# create a MembraneVector with the given mesh
pv_rbc = mir.ParticleVectors.MembraneVector("rbc", mass=1.0, mesh=mesh_rbc)
# place initial membrane
# we need a position pos and an orientation described by a quaternion q
# here we create only one membrane at the center of the domain
pos_q = [
    0.5 * domain[0],
    0.5 * domain[1],
    0.5 * domain[2],  # position
    1.0,
    0.0,
    0.0,
    0.0,
]  # quaternion
ic_rbc = mir.InitialConditions.Membrane([pos_q])
u.registerParticleVector(pv_rbc, ic_rbc)

force_scale = 0.25 * ureg.pN  # arbitrary
length_scale = np.sqrt(params.A0 / A0)
time_scale = 1 * ureg.s  # arbitrary
mass_scale = (force_scale / length_scale * time_scale**2).to(ureg.kg)

rbc_params = params.get_params(
    length_scale=length_scale,
    time_scale=time_scale,
    mass_scale=mass_scale,
    mesh=mesh_ini,
)


# now we create the internal interaction
# here we take the WLC model for shear forces and Kantor model for bending forces.
# the parameters are passed in a kwargs style
int_rbc = mir.Interactions.MembraneForces("int_rbc", **rbc_params.to_interactions())

# then we proceed as usual to make th membrane particles evolve in time
vv = mir.Integrators.VelocityVerlet("vv")
u.registerIntegrator(vv)
u.setIntegrator(vv, pv_rbc)
u.registerInteraction(int_rbc)
u.setInteraction(int_rbc, pv_rbc, pv_rbc)

# dump the mesh every 50 steps in ply format to the folder 'ply/'
u.registerPlugins(
    mir.Plugins.createDumpMesh("mesh_dump", pv_rbc, 50, "membrane_param_lucas_ply/")
)

u.run(5002, dt=dt)
