#!/usr/bin/env python

import mirheo as mir

dt = 0.001
rc = 1.0
number_density = 8.0

ranks = (1, 1, 1)
domain = (16.0, 16.0, 16.0)

u = mir.Mirheo(ranks, domain, debug_level=3, log_filename="log")

# create the particle vectors
#############################

# create MembraneVector for membranes
mesh_rbc = mir.ParticleVectors.MembraneMesh("data/rbc_mesh.off")

pv_rbc = mir.ParticleVectors.MembraneVector("rbc", mass=1.0, mesh=mesh_rbc)
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
pv_outer = mir.ParticleVectors.ParticleVector("pv_outer", mass=1.0)
ic_outer = mir.InitialConditions.Uniform(number_density)
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

prms_rbc = {
    "x0": 0.457,
    "ka_tot": 4900.0,
    "kv_tot": 7500.0,
    "ka": 5000,
    "ks": 0.0444 / 0.000906667,
    "mpow": 2.0,
    "gammaC": 52.0,
    "kBT": 0.0,
    "tot_area": 62.2242,
    "tot_volume": 26.6649,
    "kb": 44.4444,
    "theta": 6.97,
}

int_rbc = mir.Interactions.MembraneForces("int_rbc", "wlc", "Kantor", **prms_rbc)
int_dpd_oo = mir.Interactions.Pairwise(
    "dpd_oo", rc, kind="DPD", a=10.0, gamma=10.0, kBT=1.0, power=0.5
)
int_dpd_ii = mir.Interactions.Pairwise(
    "dpd_ii", rc, kind="DPD", a=10.0, gamma=20.0, kBT=1.0, power=0.5
)
int_dpd_io = mir.Interactions.Pairwise(
    "dpd_io", rc, kind="DPD", a=10.0, gamma=15.0, kBT=1.0, power=0.5
)
int_dpd_sr = mir.Interactions.Pairwise(
    "dpd_sr", rc, kind="DPD", a=0.0, gamma=15.0, kBT=1.0, power=0.5
)


u.registerInteraction(int_rbc)
u.registerInteraction(int_dpd_oo)
u.registerInteraction(int_dpd_ii)
u.registerInteraction(int_dpd_io)
u.registerInteraction(int_dpd_sr)

vv = mir.Integrators.VelocityVerlet_withConstForce("vv", force=(0.2, 0, 0))
u.registerIntegrator(vv)


center = (domain[1] * 0.5, domain[2] * 0.5)  # center in the (yz) plane
radius = 0.5 * domain[1] - rc  # radius needs to be smaller than half of the domain
# because of the frozen particles

wall = mir.Walls.Cylinder(
    "cylinder", center=center, radius=radius, axis="x", inside=True
)

u.registerWall(wall)
pv_frozen = u.makeFrozenWallParticles(
    pvName="wall",
    walls=[wall],
    interactions=[int_dpd_oo],
    integrator=vv,
    number_density=number_density,
    dt=dt,
)

u.setWall(wall, pv_outer)
u.setWall(wall, pv_inner)
u.setWall(wall, pv_rbc)

u.setInteraction(int_dpd_oo, pv_outer, pv_outer)
u.setInteraction(int_dpd_oo, pv_outer, pv_frozen)
u.setInteraction(int_dpd_ii, pv_inner, pv_inner)
u.setInteraction(int_dpd_io, pv_inner, pv_outer)

u.setInteraction(int_rbc, pv_rbc, pv_rbc)

u.setInteraction(int_dpd_sr, pv_outer, pv_rbc)
u.setInteraction(int_dpd_sr, pv_inner, pv_rbc)

# integrators
#############


u.setIntegrator(vv, pv_outer)
u.setIntegrator(vv, pv_inner)
u.setIntegrator(vv, pv_rbc)

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

u.registerPlugins(mir.Plugins.createStats("stats", every=500))

dump_every = 500
u.registerPlugins(
    mir.Plugins.createDumpParticles(
        "part_dump_inner", pv_inner, dump_every, [], "h5/membrane_slovent/inner-"
    )
)
u.registerPlugins(
    mir.Plugins.createDumpParticles(
        "part_dump_outer", pv_outer, dump_every, [], "h5/membrane_slovent/outer-"
    )
)
u.registerPlugins(
    mir.Plugins.createDumpMesh(
        "mesh_dump", pv_rbc, dump_every, "h5/membrane_slovent/ply/"
    )
)

u.run(5002, dt=dt)
