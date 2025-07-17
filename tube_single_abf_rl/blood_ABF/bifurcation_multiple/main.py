#!/usr/bin/env python

import cupy as cp
import quaternion
import mirheo as mir
import numpy as np
import sys

from cylinder.parameters import (SimulationParams,
                        load_parameters)


def get_quaternion_between_vectors(u, v):
    k_cos_theta = np.dot(u, v)
    k = np.sqrt(np.dot(u, u) * np.dot(v, v))

    if k_cos_theta == -k:
        return np.array([0, 0, 1, 0])

    q = np.array([k_cos_theta + k, *np.cross(u, v)])
    return q / np.linalg.norm(q)


def run_capillary_flow(p: 'SimulationParams',
                       *,
                       restart_directory: str=None,
                       checkpoint_directory: str=None,
                       no_visc: bool=False,
                       ABF_coords: list,
                       sdf_filename: str,
                       field_angle_deg: float,
                       ranks: tuple):

    rc = p.rc
    RA = p.RA

    max_contact_force = p.rbc_params.shear_modulus() * RA * 0.5

    dt = p.dpd_ii.get_max_dt()
    tend = 100 * 2 * np.pi / p.omega
    niters = int(tend/dt)

    domain = p.domain

    uargs = {'nranks': ranks,
             'domain': domain,
             'debug_level': 3,
             'log_filename': 'log'
    }

    if checkpoint_directory is not None:
        uargs['checkpoint_every']  = 10000
        uargs['checkpoint_folder'] = checkpoint_directory

    u = mir.Mirheo(**uargs)

    vv = mir.Integrators.VelocityVerlet('vv')
    u.registerIntegrator(vv)

    # Membranes

    mesh_rbc = mir.ParticleVectors.MembraneMesh(p.mesh_ini.vertices.tolist(),
                                                p.mesh_ref.vertices.tolist(),
                                                p.mesh_ini.faces.tolist())
    pv_rbc = mir.ParticleVectors.MembraneVector('rbc', mass=p.m, mesh=mesh_rbc)
    u.registerParticleVector(pv_rbc, mir.InitialConditions.Restart('generate_ic'))

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
    u.registerObjectBelongingChecker(rbc_checker, pv_rbc)
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

    # no viscosity in any cases: no_visc or if using Shardlow integration.
    rbc_int = mir.Interactions.MembraneForces('int_rbc', **p.rbc_params.to_interactions_zero_visc(), stress_free=True)
    u.registerInteraction(rbc_int)

    sigma = RA/7
    eps = 5 * p.rbc_params.bending_modulus()
    contact = mir.Interactions.Pairwise('contact', rc, kind='RepulsiveLJ',
                                        epsilon=eps, sigma=sigma,
                                        aware_mode='Object', max_force=max_contact_force)
    u.registerInteraction(contact)

    # Walls

    wall = mir.Walls.SDF('bifurcation', sdf_filename)
    u.registerWall(wall)

    frozen = u.makeFrozenWallParticles("frozen", walls=[wall], interactions=[dpd_oo], integrator=vv,
                                       number_density=p.nd, nsteps=int(1.0/dt), dt=dt)
    u.setWall(wall, pv_outer)

    # Set interactions between pvs

    u.registerInteraction(dpd_rbco)
    u.registerInteraction(dpd_rbci)

    u.setInteraction(dpd_oo, pv_outer, pv_outer)
    u.setInteraction(dpd_ii, pv_inner, pv_inner)
    u.setInteraction(dpd_io, pv_inner, pv_outer)

    u.setInteraction(dpd_rbco, pv_rbc, pv_outer)
    u.setInteraction(dpd_rbci, pv_rbc, pv_inner)

    u.setInteraction(dpd_oo, pv_outer, frozen)
    u.setInteraction(dpd_io, pv_inner, frozen)
    u.setInteraction(dpd_rbco, pv_rbc, frozen)

    u.setInteraction(dpd_oo, pv_outer, pv_ABF)
    u.setInteraction(dpd_io, pv_inner, pv_ABF)
    u.setInteraction(dpd_rbco, pv_rbc, pv_ABF)

    u.setInteraction(contact, pv_rbc, pv_rbc)
    u.setInteraction(contact, pv_rbc, pv_ABF)
    #u.setInteraction(rbc_int, pv_rbc, pv_rbc)

    # Integrators

    dt_rbc_el = p.rbc_params.get_max_dt_elastic(mass=p.m)
    dt_rbc_visc = p.rbc_params.get_max_dt_visc(mass=p.m)
    substeps = 5 + int(dt/dt_rbc_el)

    if no_visc:
        vv_rbc = mir.Integrators.SubStep("vv_rbc", substeps=substeps,
                                         fastForces=[rbc_int])
    else:
        vv_rbc = mir.Integrators.SubStepShardlowSweep("vv_rbc", substeps=substeps,
                                                      fastForces=rbc_int,
                                                      **p.rbc_params.to_viscous(),
                                                      nsweeps=10)
    u.registerIntegrator(vv_rbc)

    u.setIntegrator(vv_rbc, pv_rbc)
    u.setIntegrator(vv, pv_outer)
    u.setIntegrator(vv, pv_inner)
    u.setIntegrator(vv_ABF, pv_ABF)

    # bouncers

    ABF_bouncer = mir.Bouncers.Mesh("ABF_bouncer", "bounce_maxwell", kBT=0.0)
    u.registerBouncer(ABF_bouncer)

    rbc_bouncer = mir.Bouncers.Mesh("rbc_bouncer", "bounce_maxwell", kBT=0.0)
    u.registerBouncer(rbc_bouncer)

    u.setBouncer(ABF_bouncer, pv_ABF, pv_outer)
    u.setBouncer(ABF_bouncer, pv_rbc, pv_outer)
    u.setBouncer(ABF_bouncer, pv_rbc, pv_inner)

    belonging_checker = mir.BelongingCheckers.Mesh("inside_ABF_checker")
    u.registerObjectBelongingChecker(belonging_checker, pv_ABF)
    u.applyObjectBelongingChecker(belonging_checker, pv_outer, correct_every=10000, inside="none", outside="")


    # Plugins

    t_dump_every = 0.04 * 2 * np.pi / p.omega

    stats_every = int(t_dump_every/dt)
    dump_every = int(t_dump_every/dt)

    h = rc
    max_force = 100 * p.dpd_oo.a
    C = max_force/h

    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force_RBC", pv_rbc, wall,
                                                      C=C, h=h, max_force=max_force))
    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force_ABF", pv_ABF, wall,
                                                      C=C, h=h, max_force=max_force))

    u.registerPlugins(mir.Plugins.createDumpMesh("RBC_mesh_dump", pv_rbc, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createDumpMesh("ABF_mesh_dump", pv_ABF, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createStats('stats.csv', every=stats_every, filename="stats.csv"))
    u.registerPlugins(mir.Plugins.createDumpObjectStats("obj_stats", ov=pv_ABF, dump_every=dump_every,
                                                        filename="obj_stats/ABF.csv"))

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
        u.registerPlugins(mir.Plugins.createVelocityControl("vel_control", "vcont.csv", [pv_rbc, pv_inner, pv_outer, pv_ABF],
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

    if restart_directory is not None:
        u.restart(restart_directory)
    else:
        u.dumpWalls2XDMF([wall], h=(0.25,0.25,0.25), filename='h5/wall')

    u.run(niters, dt)


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Run RBCs flowing in a capillary pipe.')
    parser.add_argument('parameters', type=str, help="Parameters of the simulation.")
    parser.add_argument('--restart-dir', type=str, default=None, help="The restart directory name (no restart if not set)")
    parser.add_argument('--checkpoint-dir', type=str, default=None, help="The checkpoint directory name (no checkpoint if not set)")
    parser.add_argument('--ABF-coords', type=str, required=True, help="The coordinates of the frozen particles of the ABF.")
    parser.add_argument('--sdf', type=str, required=True, help="The sdf filename.")
    parser.add_argument('--theta-deg', type=float, required=True,
                        help="The angle between the imposed swimming direction and the axis of the pipe, in degrees.")
    parser.add_argument('--ranks', type=int, nargs=3, default=(1,1,1), help="Number of ranks along each dimension.")
    parser.add_argument('--no-visc', action='store_true', default=False, help="Disable membrane viscosity.")
    args = parser.parse_args(argv)

    p = load_parameters(args.parameters)
    ABF_coords = np.loadtxt(args.ABF_coords)

    run_capillary_flow(p,
                       no_visc=args.no_visc,
                       ABF_coords=ABF_coords,
                       sdf_filename=args.sdf,
                       field_angle_deg=args.theta_deg,
                       ranks=args.ranks,
                       restart_directory=args.restart_dir,
                       checkpoint_directory=args.checkpoint_dir)

if __name__ == '__main__':
    main(sys.argv[1:])
