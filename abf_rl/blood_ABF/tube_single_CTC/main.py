#!/usr/bin/env python

import cupy as cp
import quaternion
import mirheo as mir
import numpy as np
import sys

from cylinder.parameters import (CapilarryFlowParams,
                        load_parameters)


def get_quaternion_between_vectors(u, v):
    k_cos_theta = np.dot(u, v)
    k = np.sqrt(np.dot(u, u) * np.dot(v, v))

    if k_cos_theta == -k:
        return np.array([0, 0, 1, 0])

    q = np.array([k_cos_theta + k, *np.cross(u, v)])
    return q / np.linalg.norm(q)


def run_capillary_flow(p: 'CapilarryFlowParams',
                       *,
                       ABF_coords: list,
                       CTC_coords: list,
                       control_param: float,
                       ranks: tuple=(1,1,1),
                       substeps=10):

    rc = p.rc
    L  = p.L
    R  = p.R
    RA = p.RA

    max_contact_force = p.rbc_params.shear_modulus() * RA * 0.5

    dt = p.dpd_ii.get_max_dt() / 5
    tend = 100 * 2 * np.pi / p.omega

    domain = (L, 2*R + 4*rc, 2*R + 4*rc)

    u = mir.Mirheo(ranks, domain, debug_level=3, log_filename='log')

    vv = mir.Integrators.VelocityVerlet('vv')
    u.registerIntegrator(vv)

    vv_rig = mir.Integrators.RigidVelocityVerlet("vv_rig")
    u.registerIntegrator(vv_rig)


    # Red blood cells

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

    # CTC

    R_CTC = p.CTC_diameter/2
    mesh_CTC = mir.ParticleVectors.Mesh(p.mesh_CTC.vertices.tolist(),
                                        p.mesh_CTC.faces.tolist()) # only for visualization

    pv_CTC = mir.ParticleVectors.RigidEllipsoidVector('CTC', mass=p.m, object_size=len(CTC_coords),
                                                      semi_axes=(R_CTC, R_CTC, R_CTC),
                                                      mesh=mesh_CTC)

    ic_CTC = mir.InitialConditions.Restart('generate_ic')
    u.registerParticleVector(pv_CTC, ic_CTC)


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

    rbc_int = mir.Interactions.MembraneForces('int_rbc', **p.rbc_params.to_interactions(), stress_free=True)
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

    u.setInteraction(dpd_oo, pv_outer, pv_ABF)
    u.setInteraction(dpd_io, pv_inner, pv_ABF)
    u.setInteraction(dpd_rbco, pv_rbc, pv_ABF)

    u.setInteraction(dpd_oo, pv_outer, pv_CTC)
    u.setInteraction(dpd_io, pv_inner, pv_CTC)
    u.setInteraction(dpd_rbco, pv_rbc, pv_CTC)

    u.setInteraction(contact, pv_rbc, pv_rbc)
    u.setInteraction(contact, pv_rbc, pv_ABF)
    u.setInteraction(contact, pv_rbc, pv_CTC)
    u.setInteraction(contact, pv_CTC, pv_ABF)

    # assume only one CTC and one ABF, no contact force between them.

    #u.setInteraction(rbc_int, pv_rbc, pv_rbc)

    # Integrators

    vv_rbc = mir.Integrators.SubStep("vv_rbc", substeps=substeps, fastForces=[rbc_int])
    u.registerIntegrator(vv_rbc)

    u.setIntegrator(vv_rbc, pv_rbc)
    u.setIntegrator(vv, pv_outer)
    u.setIntegrator(vv, pv_inner)
    u.setIntegrator(vv_rig, pv_ABF)
    u.setIntegrator(vv_rig, pv_CTC)

    # bouncers

    ABF_bouncer = mir.Bouncers.Mesh("ABF_bouncer", "bounce_maxwell", kBT=0.0)
    u.registerBouncer(ABF_bouncer)

    CTC_bouncer = mir.Bouncers.Mesh("CTC_bouncer", "bounce_maxwell", kBT=0.0)
    u.registerBouncer(CTC_bouncer)

    rbc_bouncer = mir.Bouncers.Mesh("rbc_bouncer", "bounce_maxwell", kBT=0.0)
    u.registerBouncer(rbc_bouncer)

    u.setBouncer(ABF_bouncer, pv_ABF, pv_outer)
    u.setBouncer(ABF_bouncer, pv_CTC, pv_outer)
    u.setBouncer(ABF_bouncer, pv_rbc, pv_outer)
    u.setBouncer(ABF_bouncer, pv_rbc, pv_inner)

    ABF_checker = mir.BelongingCheckers.Mesh("inside_ABF_checker")
    u.registerObjectBelongingChecker(ABF_checker, pv_ABF)
    u.applyObjectBelongingChecker(ABF_checker, pv_outer, correct_every=10000, inside="none", outside="")

    CTC_checker = mir.BelongingCheckers.Mesh("inside_CTC_checker")
    u.registerObjectBelongingChecker(CTC_checker, pv_CTC)
    u.applyObjectBelongingChecker(CTC_checker, pv_outer, correct_every=10000, inside="none", outside="")


    # Plugins

    t_dump_every = 0.1 * 2 * np.pi / p.omega

    stats_every = int(t_dump_every/dt)
    dump_every = int(t_dump_every/dt)

    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force_RBC", pv_rbc, wall, C=1, h=rc, max_force=max_contact_force))

    fw = 10 * p.dpd_oo.a
    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force_ABF", pv_ABF, wall, C=fw/rc, h=rc, max_force=fw))
    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force_CTC", pv_CTC, wall, C=fw/rc, h=rc, max_force=fw))

    u.registerPlugins(mir.Plugins.createDumpMesh("RBC_mesh_dump", pv_rbc, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createDumpMesh("ABF_mesh_dump", pv_ABF, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createDumpMesh("CTC_mesh_dump", pv_CTC, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createStats('stats.csv', every=stats_every))

    u.registerPlugins(mir.Plugins.createDumpObjectStats("RBC_stats", ov=pv_rbc, dump_every=dump_every, path="obj_stats"))
    u.registerPlugins(mir.Plugins.createDumpObjectStats("ABF_stats", ov=pv_ABF, dump_every=dump_every, path="obj_stats"))
    u.registerPlugins(mir.Plugins.createDumpObjectStats("CTC_stats", ov=pv_CTC, dump_every=dump_every, path="obj_stats"))

    B_magn = p.B_magn
    m_magn = p.m_magn
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

        elif state.current_step % control_update_every == 0:
            com_extents_ABF = cp.asarray(pv_ABF.local.per_object['com_extents'])
            pos = com_extents_ABF[0,0:3].get()

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
             B_magn * np.cos(omega_t),
             B_magn * np.sin(omega_t))

        return quaternion.rotate_vectors(q, B)

    moment = (0., m_magn, 0.)

    u.registerPlugins(mir.Plugins.createExternalMagneticTorque("magnetic_torque", pv_ABF, moment, magnetic_field))


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
        print(f"omega = {omega}")
        sys.stdout.flush()

    u.dumpWalls2XDMF([wall], h=(1,1,1), filename='h5/wall')

    u.run(int(tend/dt), dt)


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Run RBCs flowing in a capillary pipe.')
    parser.add_argument('parameters', type=str, help="Parameters of the simulation.")
    parser.add_argument('--ABF-coords', type=str, required=True, help="The coordinates of the frozen particles of the ABF.")
    parser.add_argument('--CTC-coords', type=str, required=True, help="The coordinates of the frozen particles of the CTC.")
    parser.add_argument('--control-param', type=float, default=0, help="The parameter used to control the ABF trajectory. 0 to disable control.")
    args = parser.parse_args(argv)

    p = load_parameters(args.parameters)
    ABF_coords = np.loadtxt(args.ABF_coords)
    CTC_coords = np.loadtxt(args.CTC_coords)

    run_capillary_flow(p,
                       ABF_coords=ABF_coords,
                       CTC_coords=CTC_coords,
                       control_param=args.control_param,
                       substeps=10)

if __name__ == '__main__':
    main(sys.argv[1:])
