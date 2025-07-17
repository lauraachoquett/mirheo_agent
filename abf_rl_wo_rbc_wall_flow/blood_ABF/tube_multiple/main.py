#!/usr/bin/env python

import mirheo as mir
import numpy as np
import sys

from cylinder.parameters import (SimulationParams,
                        load_parameters)


def run_capillary_flow(p: 'SimulationParams',
                       *,
                       ABF_coords: list,
                       field_angle_deg: float,
                       ranks: tuple,
                       no_visc: bool):

    rc = p.rc
    L  = p.L
    R  = p.R
    RA = p.RA

    max_contact_force = p.rbc_params.shear_modulus() * RA * 0.5

    dt = p.dpd_ii.get_max_dt()
    tend = 100 * 2 * np.pi / p.omega
    niters = int(tend/dt)

    domain = (L, 2*R + 4*rc, 2*R + 4*rc)

    u = mir.Mirheo(ranks, domain, debug_level=3, log_filename='log')

    vv = mir.Integrators.VelocityVerlet('vv')
    u.registerIntegrator(vv)

    # Membranes

    mesh_rbc = mir.ParticleVectors.MembraneMesh(p.mesh_ini.vertices.tolist(),
                                                p.mesh_ref.vertices.tolist(),
                                                p.mesh_ini.faces.tolist())
    pv_RBC = mir.ParticleVectors.MembraneVector('rbc', mass=p.m, mesh=mesh_rbc)
    u.registerParticleVector(pv_RBC, mir.InitialConditions.Restart('generate_ic'))

    # ABF

    mesh_ABF = mir.ParticleVectors.Mesh(p.mesh_ABF.vertices.tolist(), p.mesh_ABF.faces.tolist())
    I = np.diag(p.mesh_ABF.moment_inertia * p.nd * p.m).tolist()
    pv_ABF = mir.ParticleVectors.RigidObjectVector('ABF', mass=p.m, inertia=I, object_size=len(ABF_coords), mesh=mesh_ABF)
    ic_ABF = mir.InitialConditions.Restart('generate_ic')

    u.registerParticleVector(pv_ABF, ic_ABF)

    vv_ABF = mir.Integrators.RigidVelocityVerlet("vv_ABF")
    u.registerIntegrator(vv_ABF)


    # Solvent

    pv_outer = mir.ParticleVectors.ParticleVector('outer', mass=p.m)
    ic_outer = mir.InitialConditions.Uniform(number_density=p.nd)
    u.registerParticleVector(pv=pv_outer, ic=ic_outer)

    rbc_checker = mir.BelongingCheckers.Mesh('rbc_checker')
    u.registerObjectBelongingChecker(rbc_checker, pv_RBC)
    pv_inner = u.applyObjectBelongingChecker(rbc_checker, pv_outer, inside='inner', correct_every=10000)

    # Interactions

    dpd_oo = mir.Interactions.Pairwise('dpd_oo', rc=rc, kind="DPD", **p.dpd_oo.to_interactions())
    dpd_ii = mir.Interactions.Pairwise('dpd_ii', rc=rc, kind="DPD", **p.dpd_ii.to_interactions())
    dpd_io = mir.Interactions.Pairwise('dpd_io', rc=rc, kind="DPD", **p.dpd_io.to_interactions())

    u.registerInteraction(dpd_oo)
    u.registerInteraction(dpd_ii)
    u.registerInteraction(dpd_io)

    dpd_rbco = mir.Interactions.Pairwise('dpd_rbco', rc=rc, kind="DPD", **p.dpd_rbco.to_interactions())
    dpd_rbci = mir.Interactions.Pairwise('dpd_rbci', rc=rc, kind="DPD", **p.dpd_rbci.to_interactions())

    # always zero viscous interactions here: either in shardlow integrator or zero visc
    rbc_int = mir.Interactions.MembraneForces('int_rbc', **p.rbc_params.to_interactions_zero_visc(), stress_free=True)
    u.registerInteraction(rbc_int)

    sigma = RA/7
    eps = 5 * p.rbc_params.bending_modulus()
    contact = mir.Interactions.Pairwise('contact', rc, kind='RepulsiveLJ',
                                        epsilon=eps, sigma=sigma, aware_mode='Object', max_force=max_contact_force)
    u.registerInteraction(contact)

    # Walls

    wall = mir.Walls.Cylinder('pipe', center=(domain[1]/2, domain[2]/2),
                              radius=R, axis='x', inside=True)
    u.registerWall(wall)
    frozen = u.makeFrozenWallParticles("frozen", walls=[wall], interactions=[dpd_oo], integrator=vv, number_density=p.nd, nsteps=int(1.0/dt), dt=dt)
    u.setWall(wall, pv_outer)

    # Set interactions between pvs

    u.registerInteraction(dpd_rbco)
    u.registerInteraction(dpd_rbci)

    u.setInteraction(dpd_oo, pv_outer, pv_outer)
    u.setInteraction(dpd_ii, pv_inner, pv_inner)
    u.setInteraction(dpd_io, pv_inner, pv_outer)

    u.setInteraction(dpd_rbco, pv_RBC, pv_outer)
    u.setInteraction(dpd_rbci, pv_RBC, pv_inner)

    u.setInteraction(dpd_oo, pv_outer, frozen)
    u.setInteraction(dpd_io, pv_inner, frozen)
    u.setInteraction(dpd_rbco, pv_RBC, frozen)
    u.setInteraction(dpd_oo,   pv_ABF, frozen)

    u.setInteraction(dpd_oo, pv_outer, pv_ABF)
    u.setInteraction(dpd_io, pv_inner, pv_ABF)
    u.setInteraction(dpd_rbco, pv_RBC, pv_ABF)
    u.setInteraction(dpd_oo,   pv_ABF, pv_ABF)

    u.setInteraction(contact, pv_RBC, pv_RBC)
    u.setInteraction(contact, pv_RBC, pv_ABF)
    u.setInteraction(contact, pv_ABF, pv_ABF)
    #u.setInteraction(rbc_int, pv_RBC, pv_RBC)

    # Integrators

    dt_rbc_el = p.rbc_params.get_max_dt_elastic(mass=p.m)
    dt_rbc_visc = p.rbc_params.get_max_dt_visc(mass=p.m)
    substeps = 5 + int(dt/dt_rbc_el)

    if no_visc:
        vv_rbc = mir.Integrators.SubStep("vv_rbc", substeps=substeps, fastForces=[rbc_int])
    else:
        vv_rbc = mir.Integrators.SubStepShardlowSweep("vv_rbc", substeps=substeps,
                                                      fastForces=rbc_int, **p.rbc_params.to_viscous(),
                                                      nsweeps=10)
    u.registerIntegrator(vv_rbc)

    u.setIntegrator(vv_rbc, pv_RBC)
    u.setIntegrator(vv, pv_outer)
    u.setIntegrator(vv, pv_inner)
    u.setIntegrator(vv_ABF, pv_ABF)

    # bouncers

    ABF_bouncer = mir.Bouncers.Mesh("ABF_bouncer", "bounce_maxwell", kBT=0.0)
    u.registerBouncer(ABF_bouncer)

    rbc_bouncer = mir.Bouncers.Mesh("rbc_bouncer", "bounce_maxwell", kBT=0.0)
    u.registerBouncer(rbc_bouncer)

    u.setBouncer(ABF_bouncer, pv_ABF, pv_outer)
    u.setBouncer(ABF_bouncer, pv_RBC, pv_outer)
    u.setBouncer(ABF_bouncer, pv_RBC, pv_inner)

    belonging_checker = mir.BelongingCheckers.Mesh("inside_ABF_checker")
    u.registerObjectBelongingChecker(belonging_checker, pv_ABF)
    u.applyObjectBelongingChecker(belonging_checker, pv_outer, correct_every=10000, inside="none", outside="")


    # Plugins

    t_dump_every = 0.04 * 2 * np.pi / p.omega

    stats_every = int(t_dump_every/dt)
    dump_every = int(t_dump_every/dt)

    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force_RBC", pv_RBC, wall, C=1, h=rc, max_force=max_contact_force))
    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force_ABF", pv_ABF, wall, C=10*p.dpd_oo.a/rc, h=rc, max_force=10*p.dpd_oo.a))

    u.registerPlugins(mir.Plugins.createDumpMesh("RBC_mesh_dump", pv_RBC, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createDumpMesh("ABF_mesh_dump", pv_ABF, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createStats('stats.csv', every=stats_every))
    u.registerPlugins(mir.Plugins.createDumpObjectStats("obj_stats_abf", ov=pv_ABF,
                                                        dump_every=dump_every, filename="stats/ABF.csv"))
    u.registerPlugins(mir.Plugins.createDumpObjectStats("obj_stats_rbc", ov=pv_RBC,
                                                        dump_every=dump_every, filename="stats/RBC.csv"))

    magn_B = p.magn_B
    omega = p.omega
    theta = field_angle_deg * np.pi / 180

    def magnetic_field(t):
        omega_t = omega * t
        B = (magn_B * np.cos(omega_t) * np.sin(theta),
             magn_B * np.cos(omega_t) * np.cos(theta),
             magn_B * np.sin(omega_t))
        return B

    u.registerPlugins(mir.Plugins.createExternalMagneticTorque("magnetic_torque", pv_ABF, p.magn_m, magnetic_field))

    mean_vel = p.mean_vel
    if mean_vel > 0:
        factor = 0.08
        Kp = 2.0 * factor
        Ki = 1.0 * factor
        Kd = 8.0 * factor
        u.registerPlugins(mir.Plugins.createVelocityControl("vel_control", "vcont.csv", [pv_RBC, pv_inner, pv_outer, pv_ABF],
                                                            low=(0,0,0), high=domain, sample_every=1, tune_every=50,
                                                            dump_every=dump_every, target_vel=(mean_vel, 0, 0),
                                                            Kp=Kp, Ki=Ki, Kd=Kd))



    if u.isMasterTask():
        print(f"tend = {tend}")
        print(f"dt = {dt}")
        print(f"substeps = {substeps}")
        print(f"omega = {omega}")
        print(f"dump every = {dump_every}")
        sys.stdout.flush()

    u.dumpWalls2XDMF([wall], h=(1,1,1), filename='h5/wall')

    u.run(niters, dt)


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Run RBCs flowing in a capillary pipe.')
    parser.add_argument('parameters', type=str, help="Parameters of the simulation.")
    parser.add_argument('--ABF-coords', type=str, required=True, help="The coordinates of the frozen particles of the ABF.")
    parser.add_argument('--theta-deg', type=float, required=True, help="The angle between the imposed swimming direction and the axis of the pipe, in degrees.")
    parser.add_argument('--no-visc', action='store_true', default=False, help="Deactivate membrane viscosity.")
    parser.add_argument('--ranks', type=int, nargs=3, default=[1,1,1], help="Number of ranks along each space dimension.")
    args = parser.parse_args(argv)

    p = load_parameters(args.parameters)
    ABF_coords = np.loadtxt(args.ABF_coords)

    run_capillary_flow(p,
                       ABF_coords=ABF_coords,
                       field_angle_deg=args.theta_deg,
                       ranks=args.ranks,
                       no_visc=args.no_visc)

if __name__ == '__main__':
    main(sys.argv[1:])
