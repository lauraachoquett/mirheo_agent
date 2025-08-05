#!/usr/bin/env python

import argparse
import mirheo as mir
from mpi4py import MPI
import numpy as np
import sys

from parameters import (load_parameters,
                        Parameters)


def run(*,
        p: Parameters,
        ABF_coords: list,
        comm: MPI.Comm,
        ranks: tuple=(1, 1, 1)):

    m  = p.mass
    nd = p.nd
    rc = p.rc
    dt = p.dpd_oo.get_max_dt() / 5

    domain = p.domain

    u = mir.Mirheo(ranks, domain,
                   debug_level=0, log_filename='log', no_splash=True,
                   comm_ptr=MPI._addressof(comm))


    # ABF

    ABF_position = np.array(domain) / 2
    mesh_ABF = mir.ParticleVectors.Mesh(p.ABF_mesh.vertices.tolist(), p.ABF_mesh.faces.tolist())
    I = np.diag(p.ABF_mesh.moment_inertia * nd * m).tolist()
    pv_ABF = mir.ParticleVectors.RigidObjectVector('ABF', mass=m, inertia=I,
                                                   object_size=len(ABF_coords),
                                                   mesh=mesh_ABF)
    ic_ABF = mir.InitialConditions.Rigid([ABF_position.tolist() + [1, 0, 0, 0]], ABF_coords)

    u.registerParticleVector(pv=pv_ABF, ic=ic_ABF)

    # Solvent

    pv_sol = mir.ParticleVectors.ParticleVector('solvent', mass=m)
    ic_sol = mir.InitialConditions.Uniform(number_density=nd)
    u.registerParticleVector(pv=pv_sol, ic=ic_sol)


    # Interactions

    dpd_oo = mir.Interactions.Pairwise('dpd_oo', rc=rc, kind="DPD", **p.dpd_oo.to_interactions())
    u.registerInteraction(dpd_oo)

    u.setInteraction(dpd_oo, pv_sol, pv_sol)
    u.setInteraction(dpd_oo, pv_ABF, pv_sol)


    # Integrators

    vv_ABF = mir.Integrators.RigidVelocityVerlet("vv_ABF")
    u.registerIntegrator(vv_ABF)
    u.setIntegrator(vv_ABF, pv_ABF)

    vv = mir.Integrators.VelocityVerlet('vv')
    u.registerIntegrator(vv)
    u.setIntegrator(vv, pv_sol)


    # Bouncers

    ABF_bouncer = mir.Bouncers.Mesh("ABF_bouncer", "bounce_maxwell", kBT=p.dpd_oo.kBT)
    u.registerBouncer(ABF_bouncer)
    u.setBouncer(ABF_bouncer, pv_ABF, pv_sol)

    belonging_checker = mir.BelongingCheckers.Mesh("inside_ABF_checker")
    u.registerObjectBelongingChecker(belonging_checker, pv_ABF)
    u.applyObjectBelongingChecker(belonging_checker, pv_sol,
                                  correct_every=10000,
                                  inside="none", outside="")


    # Plugins

    t_dump_every = 0.1 * 2 * np.pi / p.omega

    stats_every = int(t_dump_every/dt)
    dump_every  = int(t_dump_every/dt)

    u.registerPlugins(mir.Plugins.createDumpMesh("ABF_mesh_dump", pv_ABF, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createStats('stats.csv', every=stats_every))
    u.registerPlugins(mir.Plugins.createDumpObjectStats("obj_stats", ov=pv_ABF,
                                                        dump_every=dump_every,
                                                        filename="obj_stats/ABF.csv"))


    magn_B = p.magn_B
    magn_m = p.magn_m
    omega  = p.omega

    def magnetic_field(t):
        omega_t = omega * t
        return (0.0,
                magn_B * np.cos(omega_t),
                magn_B * np.sin(omega_t))

    u.registerPlugins(mir.Plugins.createExternalMagneticTorque("magnetic_torque",
                                                               pv_ABF,
                                                               magn_m,
                                                               magnetic_field))

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
    parser.add_argument('--ABF-coords', type=str, required=True,
                        help="The coordinates of the frozen particles of the ABF.")
    args = parser.parse_args(argv)

    p = load_parameters(args.parameters)
    ABF_coords = np.loadtxt(args.ABF_coords)

    comm = MPI.COMM_WORLD

    run(p=p,
        ABF_coords=ABF_coords,
        comm=comm,
        ranks=args.ranks)

if __name__ == '__main__':
    main(sys.argv[1:])
