#!/usr/bin/env python

import mirheo as mir
import numpy as np
import sys

from cylinder.parameters import (CapilarryFlowParams,
                        load_parameters)
import objplacement

def gen_com_q_pipe(*,
                   ABF_position: np.ndarray,
                   CTC_position: np.ndarray,
                   domain: tuple,
                   L: float,
                   R: float,
                   RBC_mesh,
                   ABF_mesh,
                   CTC_diameter: float,
                   scale_ini: float,
                   Ht: float,
                   seed: int=4242):
    """
    Create initial positions for RBCs in a pipe, scaled to a smaller mesh.
    """
    RBC_extents = np.ptp(RBC_mesh.vertices, axis=0) * scale_ini
    ABF_extents = np.ptp(ABF_mesh.vertices, axis=0)
    CTC_extents = np.array([CTC_diameter, CTC_diameter, CTC_diameter])

    base_geom = objplacement.Cylinder(L=domain,
                                      radius=R,
                                      center=np.array(domain)/2,
                                      axis=0)

    ABF_positions = ABF_position.reshape((-1,3))
    CTC_positions = CTC_position.reshape((-1,3))

    geom_ABF = objplacement.GeometryWithHoles(base=base_geom,
                                              hole_positions=ABF_positions,
                                              hole_extents=ABF_extents)

    geom = objplacement.GeometryWithHoles(base=geom_ABF,
                                          hole_positions=CTC_positions,
                                          hole_extents=CTC_extents)

    com_q = objplacement.generate_ic(geometry=geom,
                                     obj_volume=RBC_mesh.volume,
                                     obj_extents=RBC_extents,
                                     target_volume_fraction=Ht,
                                     seed=seed)

    return com_q


def generate_cells(p,
                   *,
                   ABF_coords: list,
                   CTC_coords: list,
                   Ht: float,
                   ABF_uncenter_factor: float,
                   ranks: tuple=(1,1,1),
                   seed: int=4242):

    rc = p.rc
    L  = p.L
    R  = p.R
    RA = p.RA
    kb = p.rbc_params.bending_modulus()

    assert ABF_uncenter_factor <= 1
    assert ABF_uncenter_factor >= 0

    max_contact_force = p.rbc_params.shear_modulus() * RA * 10
    drag = 2

    domain = (L, 2*R + 4*rc, 2*R + 4*rc)

    ABF_position = np.array([domain[0]/4, domain[1]/2, domain[2]/2])
    CTC_position = np.array([3*domain[0]/4, domain[1]/2, domain[2]/2])

    ABF_radius = p.ABF_radius
    ABF_max_radial_pos = R - ABF_radius * 1.1 # safety
    ABF_position[1] -= ABF_uncenter_factor * ABF_max_radial_pos

    scale_ini = 0.5
    com_q = gen_com_q_pipe(ABF_position=ABF_position,
                           CTC_position=CTC_position,
                           domain=domain,
                           L=L, R=R,
                           RBC_mesh=p.mesh_ini,
                           ABF_mesh=p.mesh_ABF,
                           CTC_diameter=p.CTC_diameter,
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

    contact = mir.Interactions.Pairwise('contact', rc, kind='RepulsiveLJ', epsilon=eps, sigma=sigma, aware_mode='Object', max_force=max_contact_force)
    u.registerInteraction(contact)


    p.rbc_params.gamma = 300
    rbc_int = mir.Interactions.MembraneForces('int_rbc', **p.rbc_params.to_interactions(), stress_free=True,
                                              grow_until = tend*0.5, init_length_fraction=scale_ini)
    u.registerInteraction(rbc_int)

    vv = mir.Integrators.VelocityVerlet('vv')
    u.registerIntegrator(vv)

    wall = mir.Walls.Cylinder('pipe', center=(domain[1]/2, domain[2]/2),
                              radius=R, axis='x', inside=True)
    u.registerWall(wall)

    mesh_rbc = mir.ParticleVectors.MembraneMesh(p.mesh_ini.vertices.tolist(),
                                                p.mesh_ref.vertices.tolist(),
                                                p.mesh_ini.faces.tolist())
    pv_rbc = mir.ParticleVectors.MembraneVector('rbc', mass=p.m, mesh=mesh_rbc)
    u.registerParticleVector(pv_rbc, mir.InitialConditions.Membrane(com_q, global_scale=scale_ini))

    mesh_ABF = mir.ParticleVectors.Mesh(p.mesh_ABF.vertices.tolist(),
                                        p.mesh_ABF.faces.tolist())
    I = np.diag(p.mesh_ABF.moment_inertia * p.nd * p.m).tolist()
    pv_ABF = mir.ParticleVectors.RigidObjectVector('ABF', mass=p.m, inertia=I, object_size=len(ABF_coords), mesh=mesh_ABF)
    ic_ABF = mir.InitialConditions.Rigid([ABF_position.tolist() + [1, 0, 0, 0]], ABF_coords)
    u.registerParticleVector(pv_ABF, ic_ABF)

    # model the CTC as a solid sphere
    R_CTC = p.CTC_diameter/2
    mesh_CTC = mir.ParticleVectors.Mesh(p.mesh_CTC.vertices.tolist(),
                                        p.mesh_CTC.faces.tolist()) # only for visualization

    pv_CTC = mir.ParticleVectors.RigidEllipsoidVector('CTC', mass=p.m, object_size=len(CTC_coords),
                                                      semi_axes=(R_CTC, R_CTC, R_CTC),
                                                      mesh=mesh_CTC)

    ic_CTC = mir.InitialConditions.Rigid([CTC_position.tolist() + [1, 0, 0, 0]], CTC_coords)
    u.registerParticleVector(pv_CTC, ic_CTC)

    # We fix the ABF and the CTC on purpose, no integrator for it.

    u.setInteraction(contact, pv_rbc, pv_rbc)
    u.setInteraction(contact, pv_rbc, pv_ABF)
    u.setInteraction(contact, pv_rbc, pv_CTC)

    u.setInteraction(rbc_int, pv_rbc, pv_rbc)
    u.setIntegrator(vv, pv_rbc)

    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force", pv_rbc, wall, C=1, h=rc, max_force=max_contact_force))
    u.registerPlugins(mir.Plugins.createParticleDrag("drag", pv_rbc, drag))

    if u.isMasterTask():
        print(f"domain={domain}")
        print(f"L={L}, R={R}, RA={RA}")
        sys.stdout.flush()

    u.dumpWalls2XDMF([wall], h=(1,1,1), filename='h5/wall')

    dump_every = int(5.0/dt)
    u.registerPlugins(mir.Plugins.createDumpMesh("RBC_mesh_dump", pv_rbc, dump_every, 'ply_generate_ic/'))
    u.registerPlugins(mir.Plugins.createDumpMesh("ABF_mesh_dump", pv_ABF, dump_every, 'ply_generate_ic/'))
    u.registerPlugins(mir.Plugins.createDumpMesh("CTC_mesh_dump", pv_CTC, dump_every, 'ply_generate_ic/'))
    u.registerPlugins(mir.Plugins.createStats('stats', every=dump_every))

    u.run(niters, dt)


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Generate cells with given hematocrit in pipe.')
    parser.add_argument('parameters', type=str, help="Parameters of the simulation.")
    parser.add_argument('--ABF-coords', type=str, required=True, help="The coordinates of the frozen particles of the ABF.")
    parser.add_argument('--CTC-coords', type=str, required=True, help="The coordinates of the frozen particles of the CTC.")
    parser.add_argument('--Ht', type=float, default=0.106, help="Hematocrit, in [0,1].")
    parser.add_argument('--ABF-uncenter-factor', type=float, default=0, help="How uncenter the ABF is placed, in [0,1]: 0 = centerd, 1 = at the border of the pipe")
    parser.add_argument('--seed', type=int, default=4242, help="Random seed for initial positions.")
    args = parser.parse_args(argv)

    ABF_coords = np.loadtxt(args.ABF_coords)
    CTC_coords = np.loadtxt(args.CTC_coords)

    p = load_parameters(args.parameters)

    generate_cells(p,
                   ABF_coords=ABF_coords.tolist(),
                   CTC_coords=CTC_coords.tolist(),
                   Ht=args.Ht,
                   ABF_uncenter_factor=args.ABF_uncenter_factor,
                   seed=args.seed)

if __name__ == '__main__':
    main(sys.argv[1:])
