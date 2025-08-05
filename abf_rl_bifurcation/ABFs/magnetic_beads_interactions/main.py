#!/usr/bin/env python

import argparse
import mirheo as mir
from mpi4py import MPI
import numpy as np
import sys
import trimesh

from cylinder.parameters import (load_parameters,
                        Parameters)


def run(*,
        p: Parameters,
        coords: list,
        comm: MPI.Comm,
        ranks: tuple=(1, 1, 1)):

    m  = p.mass
    nd = p.nd
    rc = p.rc
    dt = p.dpd_oo.get_max_dt()
    R = p.bead_radius

    domain = p.domain

    u = mir.Mirheo(ranks, domain,
                   debug_level=0, log_filename='log', no_splash=True,
                   comm_ptr=MPI._addressof(comm))


    # beads

    sphere = trimesh.creation.icosphere(radius=R, subdivisions=2)
    mesh_beads = mir.ParticleVectors.Mesh(vertices=sphere.vertices, faces=sphere.faces)

    pv_beads = mir.ParticleVectors.RigidEllipsoidVector('beads', mass=m,
                                                        object_size=len(coords),
                                                        semi_axes=(R, R, R),
                                                        mesh=mesh_beads)
    dr = R + 0.5 * R
    com_q = [[domain[0]/2, domain[1]/2 - dr, domain[2]/2, 1, 0, 0, 0],
             [domain[0]/2, domain[1]/2 + dr, domain[2]/2, 1, 0, 0, 0]]
    ic_beads = mir.InitialConditions.Rigid(com_q, coords)

    u.registerParticleVector(pv=pv_beads, ic=ic_beads)

    # Solvent

    pv_sol = mir.ParticleVectors.ParticleVector('solvent', mass=m)
    ic_sol = mir.InitialConditions.Uniform(number_density=nd)
    u.registerParticleVector(pv=pv_sol, ic=ic_sol)


    # Interactions

    dpd_oo = mir.Interactions.Pairwise('dpd_oo', rc=rc, kind="DPD", **p.dpd_oo.to_interactions())
    u.registerInteraction(dpd_oo)

    sigma = rc * 2**(-1/6)
    eps = 5 * p.dpd_oo.a
    max_force = 5 * p.dpd_oo.a
    contact = mir.Interactions.Pairwise('contact', rc=rc, kind="RepulsiveLJ",
                                        epsilon=eps, sigma=sigma, max_force=max_force,
                                        aware_mode="Object")
    u.registerInteraction(contact)

    u.setInteraction(dpd_oo, pv_sol, pv_sol)
    u.setInteraction(dpd_oo, pv_beads, pv_sol)
    u.setInteraction(dpd_oo, pv_beads, pv_beads)
    u.setInteraction(contact, pv_beads, pv_beads)


    # Integrators

    vv_rig = mir.Integrators.RigidVelocityVerlet("vv_rigid")
    u.registerIntegrator(vv_rig)
    u.setIntegrator(vv_rig, pv_beads)

    vv = mir.Integrators.VelocityVerlet('vv')
    u.registerIntegrator(vv)
    u.setIntegrator(vv, pv_sol)


    # Bouncers

    bouncer_beads = mir.Bouncers.Mesh("bouncer_beads", "bounce_maxwell", kBT=p.dpd_oo.kBT)
    u.registerBouncer(bouncer_beads)
    u.setBouncer(bouncer_beads, pv_beads, pv_sol)

    belonging_checker = mir.BelongingCheckers.Mesh("inside_beads_checker")
    u.registerObjectBelongingChecker(belonging_checker, pv_beads)
    u.applyObjectBelongingChecker(belonging_checker, pv_sol,
                                  correct_every=10000,
                                  inside="none", outside="")


    # Plugins

    t_dump_every = 0.05 * 2 * np.pi / p.omega

    stats_every = int(t_dump_every/dt)
    dump_every  = int(t_dump_every/dt)

    u.registerPlugins(mir.Plugins.createDumpMesh("beads_mesh_dump", pv_beads, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createStats('stats', filename='stats.csv', every=stats_every))
    u.registerPlugins(mir.Plugins.createDumpObjectStats("obj_stats", ov=pv_beads,
                                                        dump_every=dump_every,
                                                        filename="obj_stats/beads.csv"))


    magn_B = p.magn_B
    magn_m = p.magn_m
    omega  = p.omega

    def magnetic_field(t):
        omega_t = omega * t
        return (0.0,
                magn_B * np.cos(omega_t),
                magn_B * np.sin(omega_t))

    u.registerPlugins(mir.Plugins.createExternalMagneticTorque("magnetic_torque",
                                                               pv_beads,
                                                               magn_m,
                                                               magnetic_field))

    u.registerPlugins(mir.Plugins.createMagneticDipoleInteractions("magnetic_dipoles",
                                                                   pv_beads,
                                                                   magn_m,
                                                                   p.magn_permeability))

    tend = 15 * 2 * np.pi / p.omega
    niters = int(tend/dt)

    if u.isMasterTask():
        print(f"tend = {tend}")
        print(f"dt = {dt}")
        print(f"niters = {niters}")
        print(f"dump every = {dump_every}")
        sys.stdout.flush()

    u.run(niters, dt)



def main(argv):
    parser = argparse.ArgumentParser(description='Run RBCs flowing in a capillary pipe.')
    parser.add_argument('parameters', type=str, help="Parameters of the simulation.")
    parser.add_argument('--ranks', type=int, nargs=3, default=(1, 1, 1), help="Number of ranks along each direction.")
    parser.add_argument('--beads-coords', type=str, required=True,
                        help="The coordinates of the frozen particles of the bead.")
    args = parser.parse_args(argv)

    p = load_parameters(args.parameters)
    coords = np.loadtxt(args.beads_coords)

    comm = MPI.COMM_WORLD

    run(p=p,
        coords=coords,
        comm=comm,
        ranks=args.ranks)

if __name__ == '__main__':
    main(sys.argv[1:])
