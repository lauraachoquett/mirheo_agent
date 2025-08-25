#!/usr/bin/env python

import mirheo as mir
import numpy as np
import sys

from parameters import (Parameters,
                        ContactParams,
                        load_parameters)


def run(p: Parameters,
        *,
        no_visc: bool,
        ranks: tuple=(1,1,1),
        restart_directory: str=None,
        checkpoint_directory: str=None):
    """
    Arguments:
        p: the simulation parameters (see parameters.py)
        ranks: number of ranks per dimension
        restart_directory: if set, specify from which directory to restart.
        checkpoint_directory: if set, specify to which directory to dump checkpoints.
    """

    rc     = p.rc
    domain = p.domain
    RA     = p.RA

    dt_fluid = min([dpd.get_max_dt() for dpd in [p.dpd_ii, p.dpd_oo]]) / 5

    tend = 1000 * RA / p.U

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

    # Solvent

    pv_outer = mir.ParticleVectors.ParticleVector('outer', mass=p.m)
    ic_outer = mir.InitialConditions.Uniform(number_density=p.nd)
    u.registerParticleVector(pv=pv_outer, ic=ic_outer)

    rbc_checker = mir.BelongingCheckers.Mesh('rbc_checker')
    u.registerObjectBelongingChecker(rbc_checker, pv_rbc)
    pv_inner = u.applyObjectBelongingChecker(rbc_checker, pv_outer,
                                             inside='inner', correct_every=10000)

    rbc_bouncer = mir.Bouncers.Mesh("rbc_bouncer", "bounce_maxwell", kBT=p.kBT)
    u.registerBouncer(rbc_bouncer)

    # Interactions

    dpd_oo = mir.Interactions.Pairwise('dpd_oo', rc=rc, kind="DPD", **p.dpd_oo.to_interactions())
    dpd_ii = mir.Interactions.Pairwise('dpd_ii', rc=rc, kind="DPD", **p.dpd_ii.to_interactions())
    dpd_io = mir.Interactions.Pairwise('dpd_io', rc=rc, kind="DPD", **p.dpd_io.to_interactions())

    u.registerInteraction(dpd_oo)
    u.registerInteraction(dpd_ii)
    u.registerInteraction(dpd_io)

    dpd_rbco = mir.Interactions.Pairwise('dpd_rbco', rc=rc, kind="DPD",
                                         **p.dpd_rbco.to_interactions())
    dpd_rbci = mir.Interactions.Pairwise('dpd_rbci', rc=rc, kind="DPD",
                                         **p.dpd_rbci.to_interactions())

    # zero viscosity here, we will either use the Shardlow integrator or not viscosity at all.
    rbc_int = mir.Interactions.MembraneForces('int_rbc',
                                              **p.rbc_params.to_interactions_zero_visc(),
                                              stress_free=True)
    u.registerInteraction(rbc_int)

    cp = p.contact_rbcs
    sigma             = cp.sigma
    eps               = cp.eps
    max_contact_force = cp.max_contact_force

    contact = mir.Interactions.Pairwise('contact', rc, kind='RepulsiveLJ',
                                        epsilon=eps, sigma=sigma, aware_mode='Object',
                                        max_force=max_contact_force)
    u.registerInteraction(contact)

    # Walls

    wall = mir.Walls.SDF('walls', p.sdf_filename)
    u.registerWall(wall)
    frozen = u.makeFrozenWallParticles("frozen", walls=[wall], interactions=[dpd_oo],
                                       integrator=vv, number_density=p.nd,
                                       nsteps=1000, dt=p.dpd_oo.get_max_dt())
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

    u.setInteraction(contact, pv_rbc, pv_rbc)

    # set bouncers
    u.setBouncer(rbc_bouncer, pv_rbc, pv_inner)
    u.setBouncer(rbc_bouncer, pv_rbc, pv_outer)

    # Integrators

    u.setIntegrator(vv, pv_outer)
    u.setIntegrator(vv, pv_inner)

    dt_rbc_el = p.rbc_params.get_max_dt_elastic(mass=p.m)
    dt_rbc_visc = p.rbc_params.get_max_dt_visc(mass=p.m)
    substeps = 5 + int(dt_fluid / dt_rbc_el)
    dt = dt_fluid

    if no_visc:
        vv_rbc = mir.Integrators.SubStep("vv_rbc", substeps=substeps,
                                         fastForces=[rbc_int])
    else:
        nsweeps = 10
        vv_rbc = mir.Integrators.SubStepShardlowSweep("vv_rbc", substeps=substeps,
                                                      fastForces=rbc_int,
                                                      **p.rbc_params.to_viscous(),
                                                      nsweeps=nsweeps)
    u.registerIntegrator(vv_rbc)
    u.setIntegrator(vv_rbc, pv_rbc)


    # Plugins

    t_dump_every = 0.5 * RA / p.U
    dump_every = int(t_dump_every/dt)
    stats_every = dump_every

    h = p.wall_repulsion_length
    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force", pv_rbc, wall,
                                                      C=max_contact_force/h, h=h,
                                                      max_force=max_contact_force))


    f = p.forces_filename
    h = (rc, rc, rc)
    u.registerPlugins(mir.Plugins.createAddForceField("body_force_solvent", pv_outer, f, h))
    u.registerPlugins(mir.Plugins.createAddForceField("body_force_hemoglobine", pv_inner, f, h))
    #u.registerPlugins(mir.Plugins.createAddForceField("body_force_membranes", pv_rbc, f, h))

    u.registerPlugins(mir.Plugins.createDumpMesh("mesh_dump", pv_rbc, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createStats("stats", every=stats_every,
                                              filename='stats.csv',
                                              pvs=[pv_inner, pv_outer, pv_rbc]))
    u.registerPlugins(mir.Plugins.createDumpObjectStats("rbc_stats", ov=pv_rbc,
                                                        dump_every=stats_every,
                                                        filename="obj_stats/rbc.csv"))

    u.registerPlugins(mir.Plugins.createDumpAverage(name="field_dump",
                                                    pvs=[pv_outer, pv_inner],
                                                    sample_every=2,
                                                    dump_every=dump_every,
                                                    bin_size=(1.0, 1.0, 1.0),
                                                    channels=["velocities"],
                                                    path='h5/field'))



    if u.isMasterTask():
        print(f"tend = {tend}")
        print(f"dt = {dt}")
        print(f"substeps = {substeps}")
        if not no_visc:
            print(f"nsweeps = {nsweeps}")
        sys.stdout.flush()

    if restart_directory is not None:
        u.restart(restart_directory)
    else:
        u.dumpWalls2XDMF([wall], h=(1,1,1), filename='h5/wall')

    u.run(int(tend/dt), dt)


def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--params',         type=str, default="parameters.pkl", help="Parameters file.")
    parser.add_argument('--ranks',          type=int, nargs=3, default=(1, 1, 1), help="Ranks in each dimension.")
    parser.add_argument('--restart-dir',    type=str, default=None, help="The restart directory name (no restart if not set)")
    parser.add_argument('--checkpoint-dir', type=str, default=None, help="The checkpoint directory name (no checkpoint if not set)")
    parser.add_argument('--no-visc', action='store_true', default=False, help="Disable viscosity")
    args = parser.parse_args(argv)

    p = load_parameters(args.params)

    run(p, ranks=tuple(args.ranks),
        no_visc=args.no_visc,
        restart_directory=args.restart_dir,
        checkpoint_directory=args.checkpoint_dir)


if __name__ == '__main__':
    main(sys.argv[1:])
