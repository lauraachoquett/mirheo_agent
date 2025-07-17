#!/usr/bin/env python

try:
    import cupy as cp
except Exception as e:
    print("cupy is not installed. This will cause an error if control is set.")
    cp = None

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
                       abf_coords: list,
                       control_param: float,
                       ranks: tuple=(1,1,1),
                       restart_directory: str=None,
                       checkpoint_directory: str=None):

    rc = p.rc
    L  = p.L
    R  = p.R
    RA = p.RA

    max_contact_force = p.rbc_params.shear_modulus() * RA * 0.5

    dt = p.dpd_ii.get_max_dt()
    tend = 100 * 2 * np.pi / p.omega
    niters = int(tend/dt)

    domain = (L, 2*R + 4*rc, 2*R + 4*rc)

    uargs = {
        'nranks': ranks,
        'domain': domain,
        'debug_level': 0,
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

    mesh_abf = mir.ParticleVectors.Mesh(p.mesh_abf.vertices.tolist(), p.mesh_abf.faces.tolist())
    I = np.diag(p.mesh_abf.moment_inertia * p.nd * p.m).tolist()
    pv_abf = mir.ParticleVectors.RigidObjectVector('abf', mass=p.m, inertia=I, object_size=len(abf_coords), mesh=mesh_abf)
    ic_abf = mir.InitialConditions.Restart('generate_ic')

    u.registerParticleVector(pv_abf, ic_abf)

    vv_abf = mir.Integrators.RigidVelocityVerlet("vv_abf")
    u.registerIntegrator(vv_abf)


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

    rbc_int = mir.Interactions.MembraneForces('int_rbc', **p.rbc_params.to_interactions_zero_visc(), stress_free=True)
    u.registerInteraction(rbc_int)

    sigma = RA/7
    eps = 5 * p.rbc_params.bending_modulus()
    contact = mir.Interactions.Pairwise('contact', rc, kind='RepulsiveLJ', epsilon=eps, sigma=sigma, aware_mode='Object', max_force=max_contact_force)
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

    u.setInteraction(dpd_rbco, pv_rbc, pv_outer)
    u.setInteraction(dpd_rbci, pv_rbc, pv_inner)

    u.setInteraction(dpd_oo, pv_outer, frozen)
    u.setInteraction(dpd_io, pv_inner, frozen)
    u.setInteraction(dpd_rbco, pv_rbc, frozen)

    u.setInteraction(dpd_oo, pv_outer, pv_abf)
    u.setInteraction(dpd_io, pv_inner, pv_abf)
    u.setInteraction(dpd_rbco, pv_rbc, pv_abf)

    u.setInteraction(contact, pv_rbc, pv_rbc)
    u.setInteraction(contact, pv_rbc, pv_abf)
    u.setInteraction(rbc_int, pv_rbc, pv_rbc)

    # Integrators

    dt_rbc_el = p.rbc_params.get_max_dt_elastic(mass=p.m)
    dt_rbc_visc = p.rbc_params.get_max_dt_visc(mass=p.m)
    substeps = 5 + int(dt/dt_rbc_el)

    vv_rbc = mir.Integrators.SubStepShardlowSweep("vv_rbc", substeps=substeps,
                                                  fastForces=rbc_int, **p.rbc_params.to_viscous(),
                                                  nsweeps=10)
    u.registerIntegrator(vv_rbc)

    u.setIntegrator(vv_rbc, pv_rbc)
    u.setIntegrator(vv, pv_outer)
    u.setIntegrator(vv, pv_inner)
    u.setIntegrator(vv_abf, pv_abf)

    # bouncers

    abf_bouncer = mir.Bouncers.Mesh("abf_bouncer", "bounce_maxwell", kBT=0.0)
    u.registerBouncer(abf_bouncer)

    rbc_bouncer = mir.Bouncers.Mesh("rbc_bouncer", "bounce_maxwell", kBT=0.0)
    u.registerBouncer(rbc_bouncer)

    u.setBouncer(abf_bouncer, pv_abf, pv_outer)
    u.setBouncer(rbc_bouncer, pv_rbc, pv_outer)
    u.setBouncer(rbc_bouncer, pv_rbc, pv_inner)

    belonging_checker = mir.BelongingCheckers.Mesh("inside_abf_checker")
    u.registerObjectBelongingChecker(belonging_checker, pv_abf)
    u.applyObjectBelongingChecker(belonging_checker, pv_outer, correct_every=10000, inside="none", outside="")


    # Plugins

    t_dump_every = 0.04 * 2 * np.pi / p.omega

    stats_every = int(t_dump_every/dt)
    dump_every = int(t_dump_every/dt)

    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force_rbc", pv_rbc, wall, C=1, h=rc, max_force=max_contact_force))
    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force_abf", pv_abf, wall, C=10*p.dpd_oo.a/rc, h=rc, max_force=10*p.dpd_oo.a))

    u.registerPlugins(mir.Plugins.createDumpMesh("rbc_mesh_dump", pv_rbc, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createDumpMesh("abf_mesh_dump", pv_abf, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createStats('stats.csv', every=stats_every))
    u.registerPlugins(mir.Plugins.createDumpObjectStats("obj_stats", ov=pv_abf, dump_every=dump_every, filename="obj_stats/abf.csv"))

    magn_B = p.magn_B
    omega = p.omega


    state = u.getState()
    t_control_update_every = 2 * np.pi / p.omega
    control_update_every = int(t_control_update_every / dt)

    def magnetic_field(t):
        omega_t = omega * t
        ex = (1, 0, 0)

        global B_direction

        if control_param == 0:
            B_direction = (1, 0, 0)

        elif not "B_direction" in locals() or state.current_step % control_update_every == 0:
            com_extents_abf = cp.asarray(pv_abf.local.per_object['com_extents'])
            pos = com_extents_abf[0,0:3].get()

            distance_to_centerline = (pos[1]**2 + pos[2]**2)**0.5
            theta = - np.arctan(control_param * distance_to_centerline / R)
            r = [pos[1] / (1e-6 + distance_to_centerline),
                 pos[2] / (1e-6 + distance_to_centerline)]

            B_direction = (np.cos(theta),
                           r[0] * np.sin(theta),
                           r[1] * np.sin(theta))
            #print(B_direction)

        q = get_quaternion_between_vectors(ex, B_direction)
        q = quaternion.from_float_array(q)

        B = (0.0,
             magn_B * np.cos(omega_t),
             magn_B * np.sin(omega_t))

        return quaternion.rotate_vectors(q, B)

    u.registerPlugins(mir.Plugins.createExternalMagneticTorque("magnetic_torque", pv_abf, p.magn_m, magnetic_field))

    mean_vel = p.mean_vel
    if mean_vel > 0:
        factor = 0.08
        Kp = 2.0 * factor
        Ki = 1.0 * factor
        Kd = 8.0 * factor
        u.registerPlugins(mir.Plugins.createVelocityControl("vel_control", "vcont.csv", [pv_rbc, pv_inner, pv_outer, pv_abf],
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
        u.dumpWalls2XDMF([wall], h=(1,1,1), filename='h5/wall')

    u.run(niters, dt)


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Run rbcs flowing in a capillary pipe.')
    parser.add_argument('parameters', type=str, help="Parameters of the simulation.")
    parser.add_argument('--abf-coords', type=str, required=True, help="The coordinates of the frozen particles of the abf.")
    parser.add_argument('--control-param', type=float, default=0, help="The parameter used to control the abf trajectory. 0 to disable control.")
    parser.add_argument('--restart-dir', type=str, default=None, help="The restart directory name (no restart if not set)")
    parser.add_argument('--checkpoint-dir', type=str, default=None, help="The checkpoint directory name (no checkpoint if not set)")
    args = parser.parse_args(argv)

    p = load_parameters(args.parameters)
    abf_coords = np.loadtxt(args.abf_coords)

    run_capillary_flow(p,
                       abf_coords=abf_coords,
                       control_param=args.control_param,
                       restart_directory=args.restart_dir,
                       checkpoint_directory=args.checkpoint_dir)

if __name__ == '__main__':
    main(sys.argv[1:])
