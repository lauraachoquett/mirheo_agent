#!/usr/bin/env python
from copy import deepcopy
from cylinder.parameters import set_parameters
import mirheo as mir
from dpdprops import (
    JuelicherLimRBCDefaultParams,
    create_dpd_params_from_props,
    create_fsi_dpd_params,
    load_stress_free_mesh,
    load_equilibrium_mesh,
    sound_speed,
    create_dpd_params_from_Re_Ma,
)
import pint
import trimesh
import numpy as np

params = set_parameters()
rbc_params = params["rbc_params"]
params_outer = params["params_outer"]
params_inner = params["params_inner"]
params_inner_outer = params["params_inner_outer"]
inner_rbc_params = params["inner_rbc_params"]
outer_rbc_params = params["outer_rbc_params"]
contact_rbcs = params["contact_rbcs"]
dt0 = params["dt0"]
mesh_ini = params["mesh_ini"]
mesh_ref = params["mesh_ref"]
nd = params["nd"]
mass = params["mass"]
rc = params["rc"]
ranks = (2, 1, 1)
domain = (16.0, 16.0, 16.0)

u = mir.Mirheo(ranks, domain, debug_level=3, log_filename="log")
print("RBC parameters :", rbc_params)
# create the particle vectors
#############################

# create MembraneVector for membranes
mesh_rbc = mir.ParticleVectors.MembraneMesh(
    mesh_ini.vertices.tolist(), mesh_ref.vertices.tolist(), mesh_ini.faces.tolist()
)


pv_rbc = mir.ParticleVectors.MembraneVector("rbc", mass=mass, mesh=mesh_rbc)
pos_q = [
    0.3 * domain[0],
    0.5 * domain[1],
    0.5 * domain[2],  # position
    1.0,
    0.0,
    0.0,
    0.0,
]  # quaternion
ic_rbc = mir.InitialConditions.Membrane([pos_q])
u.registerParticleVector(pv_rbc, ic_rbc)

# create particleVector for outer solvent
pv_outer = mir.ParticleVectors.ParticleVector("pv_outer", mass=mass)
ic_outer = mir.InitialConditions.Uniform(nd)
u.registerParticleVector(pv_outer, ic_outer)

# To create the inner solvent, we split the outer solvent (which originally occupies
# the whole domain) into outer and inner solvent
# This is done thanks to the belonging checker:
inner_checker = mir.BelongingCheckers.Mesh("inner_solvent_checker")

# the checker needs to be registered, as any other object; it is associated to a given object vector
u.registerObjectBelongingChecker(inner_checker, pv_rbc)

# we can now apply the checker to create the inner PV
pv_inner = u.applyObjectBelongingChecker(
    inner_checker, pv_outer, correct_every=0, inside="pv_inner"
)


# interactions
##############


int_rbc = mir.Interactions.MembraneForces(
    "int_rbc", **rbc_params.to_interactions_zero_visc(), stress_free=True
)


int_dpd_oo = mir.Interactions.Pairwise(
    "dpd_oo", rc, kind="DPD", **params_outer.to_interactions()
)
int_dpd_ii = mir.Interactions.Pairwise(
    "dpd_ii", rc, kind="DPD", **params_inner.to_interactions()
)
int_dpd_io = mir.Interactions.Pairwise(
    "dpd_io",
    rc,
    kind="DPD",
    a=params_inner.a,
    gamma=0,
    kBT=params_inner.kBT,
    power=params_inner.kpow,
)
int_dpd_mi = mir.Interactions.Pairwise(
    "dpd_mi", rc, kind="DPD", **inner_rbc_params.to_interactions()
)
int_dpd_mo = mir.Interactions.Pairwise(
    "dpd_mo", rc, kind="DPD", **outer_rbc_params.to_interactions()
)
# zero viscosity here, we will either use the Shardlow integrator or not viscosity at all.


u.registerInteraction(int_rbc)
u.registerInteraction(int_dpd_oo)
u.registerInteraction(int_dpd_ii)
u.registerInteraction(int_dpd_io)
u.registerInteraction(int_dpd_mi)
u.registerInteraction(int_dpd_mo)

vv = mir.Integrators.VelocityVerlet("vv")
u.registerIntegrator(vv)


# center = (domain[1]*0.5, domain[2]*0.5) # center in the (yz) plane
# radius = 0.5 * domain[1] - rc           # radius needs to be smaller than half of the domain
# because of the frozen particles

# wall = mir.Walls.Cylinder("cylinder", center=center, radius=radius, axis="x", inside=True)

# u.registerWall(wall)
# pv_frozen = u.makeFrozenWallParticles(pvName="wall", walls=[wall], interactions=[int_dpd_oo], integrator=vv, number_density=nd, dt=dt)

# u.setWall(wall, pv_outer)
# u.setWall(wall, pv_inner)
# u.setWall(wall, pv_rbc)

u.setInteraction(int_dpd_oo, pv_outer, pv_outer)
# u.setInteraction(int_dpd_oo, pv_outer, pv_frozen)
u.setInteraction(int_dpd_ii, pv_inner, pv_inner)
u.setInteraction(int_dpd_io, pv_inner, pv_outer)

u.setInteraction(int_rbc, pv_rbc, pv_rbc)

u.setInteraction(int_dpd_mo, pv_outer, pv_rbc)
u.setInteraction(int_dpd_mi, pv_inner, pv_rbc)

# integrators
#############


u.setIntegrator(vv, pv_outer)
u.setIntegrator(vv, pv_inner)

dt_fluid = min(
    params_outer.get_max_dt(),
    params_inner.get_max_dt(),
    inner_rbc_params.get_max_dt(),
    outer_rbc_params.get_max_dt(),
)

dt_rbc_el = rbc_params.get_max_dt_elastic(mass=mass)
dt_rbc_visc = rbc_params.get_max_dt_visc(mass=mass)
substeps = 5 + int(dt_fluid / dt_rbc_el)
dt = dt_fluid

vv_rbc = mir.Integrators.SubStep("vv_rbc", substeps=substeps, fastForces=[int_rbc])
u.registerIntegrator(vv_rbc)

u.setIntegrator(vv_rbc, pv_rbc)

# Bouncers
##########

# The solvent must not go through the membrane
# we can enforce this by setting a bouncer
bouncer = mir.Bouncers.Mesh("membrane_bounce", "bounce_maxwell", kBT=0.5)

# we register the bouncer object as any other object
u.registerBouncer(bouncer)

# now we can set what PVs bounce on what OV:
u.setBouncer(bouncer, pv_rbc, pv_outer)
u.setBouncer(bouncer, pv_rbc, pv_inner)

# plugins
#########

u.registerPlugins(mir.Plugins.createStats("stats", every=5000))

dump_every = 5000
u.registerPlugins(
    mir.Plugins.createDumpParticles(
        "part_dump_inner",
        pv_inner,
        dump_every,
        [],
        "h5/membrane_slovent_lucas_param/inner-",
    )
)
u.registerPlugins(
    mir.Plugins.createDumpParticles(
        "part_dump_outer",
        pv_outer,
        dump_every,
        [],
        "h5/membrane_slovent_lucas_param/outer-",
    )
)
u.registerPlugins(
    mir.Plugins.createDumpMesh(
        "mesh_dump", pv_rbc, dump_every, "h5/membrane_slovent_lucas_param/ply/"
    )
)

u.run(1000002, dt=dt)
