#!/usr/bin/env python

import mirheo as mir
import numpy as np
import sys

from cylinder.parameters import (SimulationParams,
                        load_parameters)
import objplacement

def gen_com_q_pipe(*,
                   ABF_position: np.ndarray,
                   domain: tuple,
                   RBC_mesh,
                   ABF_mesh,
                   scale_ini: float,
                   Ht: float,
                   seed: int=4242):
    """
    Create initial positions for RBCs in a pipe, scaled to a smaller mesh.
    """
    RBC_extents = np.ptp(RBC_mesh.vertices, axis=0) * scale_ini
    ABF_extents = np.ptp(ABF_mesh.vertices, axis=0)

    base_geom = objplacement.DomainGeometry(L=domain)

    ABF_positions = ABF_position.reshape((-1,3))

    geom = objplacement.GeometryWithHoles(base=base_geom,
                                          hole_positions=ABF_positions,
                                          hole_extents=ABF_extents)

    com_q = objplacement.generate_ic(geometry=geom,
                                     obj_volume=RBC_mesh.volume,
                                     obj_extents=RBC_extents,
                                     target_volume_fraction=Ht,
                                     seed=seed)

    return com_q


def generate_cells(p,
                   *,
                   ABF_coords: list,
                   sdf_filename: str,
                   Ht: float,
                   ranks: tuple,
                   seed: int):

    rc = p.rc
    R  = p.R
    RA = p.RA
    kb = p.rbc_params.bending_modulus()

    max_contact_force = p.rbc_params.shear_modulus() * RA * 10
    drag = 2

    domain = p.domain

    ABF_radius = p.ABF_radius
    dx = domain[1] / (p.num_bif_x-1)
    ABF_position = np.array([0, (p.num_bif_x//2) * dx, domain[2]/2])

    scale_ini = 0.5
    com_q = gen_com_q_pipe(ABF_position=ABF_position,
                           domain=domain,
                           RBC_mesh=p.mesh_ini,
                           ABF_mesh=p.mesh_ABF,
                           scale_ini=scale_ini,
                           Ht=Ht, seed=seed)

    T = np.sqrt(p.m * RA**2 / kb)
    tend = 10 * T
    dt = 0.001

    niters = int(tend/dt)

    checkpoint_every = niters-5

    u = mir.Mirheo(ranks, domain, debug_level=3, log_filename='generate_ic',
                   checkpoint_folder="generate_ic/", checkpoint_every=checkpoint_every)


    sigma = RA/10
    eps = 3 * kb

    contact = mir.Interactions.Pairwise('contact', rc, kind='RepulsiveLJ',
                                        epsilon=eps, sigma=sigma,
                                        aware_mode='Object', max_force=max_contact_force)
    u.registerInteraction(contact)


    p.rbc_params.gamma = 300
    rbc_int = mir.Interactions.MembraneForces('int_rbc', **p.rbc_params.to_interactions(), stress_free=True,
                                              grow_until = tend*0.5, init_length_fraction=scale_ini)
    u.registerInteraction(rbc_int)

    vv = mir.Integrators.VelocityVerlet('vv')
    u.registerIntegrator(vv)

    wall = mir.Walls.SDF('bifurcation', sdf_filename)
    u.registerWall(wall)

    mesh_rbc = mir.ParticleVectors.MembraneMesh(p.mesh_ini.vertices.tolist(),
                                                p.mesh_ref.vertices.tolist(),
                                                p.mesh_ini.faces.tolist())
    pv_rbc = mir.ParticleVectors.MembraneVector('rbc', mass=p.m, mesh=mesh_rbc)
    u.registerParticleVector(pv_rbc, mir.InitialConditions.Membrane(com_q, global_scale=scale_ini))

    mesh_ABF = mir.ParticleVectors.Mesh(p.mesh_ABF.vertices.tolist(), p.mesh_ABF.faces.tolist())
    I = np.diag(p.mesh_ABF.moment_inertia * p.nd * p.m).tolist()
    pv_ABF = mir.ParticleVectors.RigidObjectVector('ABF', mass=p.m, inertia=I, object_size=len(ABF_coords), mesh=mesh_ABF)
    ic_ABF = mir.InitialConditions.Rigid([ABF_position.tolist() + [1, 0, 0, 0]], ABF_coords)
    # We fix th ABF on purpose, no integrator for it.
    #vv_ABF = mir.Integrators.RigidVelocityVerlet("vv_ABF")
    u.registerParticleVector(pv_ABF, ic_ABF)

    u.setInteraction(contact, pv_rbc, pv_rbc)
    u.setInteraction(contact, pv_rbc, pv_ABF)

    u.setInteraction(rbc_int, pv_rbc, pv_rbc)
    u.setIntegrator(vv, pv_rbc)

    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force", pv_rbc, wall, C=1, h=rc, max_force=max_contact_force))
    u.registerPlugins(mir.Plugins.createParticleDrag("drag", pv_rbc, drag))

    if u.isMasterTask():
        print(f"domain={domain}")
        print(f"RA={RA}")
        sys.stdout.flush()

    u.dumpWalls2XDMF([wall], h=(1,1,1), filename='h5/wall')

    dump_every = int(niters//5)
    u.registerPlugins(mir.Plugins.createDumpMesh("RBC_mesh_dump", pv_rbc, dump_every, 'ply_generate_ic/'))
    u.registerPlugins(mir.Plugins.createDumpMesh("ABF_mesh_dump", pv_ABF, dump_every, 'ply_generate_ic/'))
    u.registerPlugins(mir.Plugins.createStats('stats', every=dump_every))

    u.run(niters, dt)


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Generate cells with given hematocrit in pipe.')
    parser.add_argument('parameters', type=str, help="Parameters of the simulation.")
    parser.add_argument('--ABF-coords', type=str, required=True, help="The coordinates of the frozen particles of the ABF.")
    parser.add_argument('--Ht', type=float, default=0.106, help="Hematocrit, in [0,1].")
    parser.add_argument('--sdf', type=str, required=True, help="sdf file.")
    parser.add_argument('--seed', type=int, default=4242, help="Random seed for initial positions.")
    parser.add_argument('--ranks', type=int, nargs=3, default=(1,1,1), help="Number of ranks along each dimension.")
    args = parser.parse_args(argv)

    ABF_coords = np.loadtxt(args.ABF_coords)

    p = load_parameters(args.parameters)

    generate_cells(p,
                   ABF_coords=ABF_coords.tolist(),
                   Ht=args.Ht,
                   sdf_filename=args.sdf,
                   ranks=args.ranks,
                   seed=args.seed)

if __name__ == '__main__':
    main(sys.argv[1:])
