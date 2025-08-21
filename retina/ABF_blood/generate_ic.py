#!/usr/bin/env python

import mirheo as mir
import numpy as np
import sys

from parameters import (Parameters,
                        ContactParams,
                        load_parameters)
import objplacement

def gen_com_q(*,
              abf_position: np.ndarray,
              domain: tuple,
              RBC_mesh,
              abf_mesh,
              scale_ini: float,
              Ht: float,
              seed: int):

    RBC_extents = np.ptp(RBC_mesh.vertices, axis=0) * scale_ini
    abf_extents = np.ptp(abf_mesh.vertices, axis=0)
    RBC_radius = np.max(RBC_extents) / 2

    base_geom = objplacement.DomainGeometry(L=domain)

    abf_positions = abf_position.reshape((-1,3))

    geom = objplacement.GeometryWithHoles(base=base_geom,
                                          hole_positions=abf_positions,
                                          hole_extents=abf_extents)

    com_q = objplacement.generate_ic(geometry=geom,
                                     obj_volume=RBC_mesh.volume,
                                     obj_radius=RBC_radius,
                                     target_volume_fraction=Ht,
                                     seed=seed)

    return com_q


def generate_cells(p,
                   *,
                   abf_coords: list,
                   scale_ini: float,
                   Ht: float,
                   seed: int,
                   ranks: tuple=(1,1,1),
                   position_abf : list):

    m      = p.m
    rc     = p.rc
    RA     = p.RA
    domain = p.domain

    drag = 4

    abf_position = np.array(position_abf)

    com_q = gen_com_q(abf_position=abf_position,
                      domain=domain,
                      RBC_mesh=p.mesh_ini,
                      abf_mesh=p.mesh_abf,
                      scale_ini=scale_ini,
                      Ht=Ht, seed=seed)

    T = np.sqrt(m * RA**2  / p.rbc_params.bending_modulus())
    tend = 50 * T
    dt = 0.001

    niters = int(tend/dt)
    print("Number of iterations :",niters)
    checkpoint_every = niters-5

    u = mir.Mirheo(ranks, domain, debug_level=3, log_filename='generate_ic',
                   checkpoint_folder="generate_ic/", checkpoint_every=checkpoint_every,
                   max_obj_half_length=1.5*RA)

    cp = p.contact_rbcs
    sigma             = cp.sigma
    eps               = cp.eps
    max_contact_force = cp.max_contact_force

    contact = mir.Interactions.Pairwise('contact', rc, kind='RepulsiveLJ',
                                        epsilon=eps, sigma=sigma, aware_mode='Object',
                                        max_force=max_contact_force)
    u.registerInteraction(contact)

    p.rbc_params.gamma = 300
    rbc_int = mir.Interactions.MembraneForces('int_rbc', **p.rbc_params.to_interactions(), stress_free=True,
                                              grow_until = tend*0.5, init_length_fraction=scale_ini)
    u.registerInteraction(rbc_int)

    vv = mir.Integrators.VelocityVerlet('vv')
    u.registerIntegrator(vv)

    wall = mir.Walls.SDF('walls', p.sdf_filename)
    u.registerWall(wall)

    mesh_rbc = mir.ParticleVectors.MembraneMesh(p.mesh_ini.vertices.tolist(),
                                                p.mesh_ref.vertices.tolist(),
                                                p.mesh_ini.faces.tolist())
    pv_rbc = mir.ParticleVectors.MembraneVector('rbc', mass=m, mesh=mesh_rbc)
    u.registerParticleVector(pv_rbc, mir.InitialConditions.Membrane(com_q, global_scale=scale_ini))

    mesh_abf = mir.ParticleVectors.Mesh(p.mesh_abf.vertices.tolist(), p.mesh_abf.faces.tolist())
    I = np.diag(p.mesh_abf.moment_inertia * p.nd * p.m).tolist()
    pv_abf = mir.ParticleVectors.RigidObjectVector('abf', mass=p.m, inertia=I, object_size=len(abf_coords), mesh=mesh_abf)
    q = [np.cos(np.pi/4), 0.0, 0.0, np.sin(np.pi/4)]
    ic_abf = mir.InitialConditions.Rigid([abf_position.tolist() + q], abf_coords)
    # We fix th abf on purpose, no integrator for it.
    #vv_abf = mir.Integrators.RigidVelocityVerlet("vv_abf")
    u.registerParticleVector(pv_abf, ic_abf)

    u.setInteraction(contact, pv_rbc, pv_rbc)
    u.setInteraction(contact, pv_abf, pv_rbc)

    u.setInteraction(rbc_int, pv_rbc, pv_rbc)
    u.setIntegrator(vv, pv_rbc)

    h = 0.1 * rc
    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force", pv_rbc, wall,
                                                      C=max_contact_force/h, h=h,
                                                      max_force=max_contact_force))

    u.registerPlugins(mir.Plugins.createParticleDrag("drag", pv_rbc, drag))

    # u.registerPlugins(mir.Plugins.createDumpMesh("RBC_mesh_dump", pv_rbc, checkpoint_every-100, 'ply_generate_ic/'))
    # u.registerPlugins(mir.Plugins.createDumpMesh("ABF_mesh_dump", pv_abf, checkpoint_every-100, 'ply_generate_ic/'))
    u.registerPlugins(mir.Plugins.createStats('stats', every=checkpoint_every))

    if u.isMasterTask():
        print(f"domain={domain}")
        print(f"RA={RA}")
        sys.stdout.flush()

    u.dumpWalls2XDMF([wall], h=(1,1,1), filename='h5/wall')

    u.run(niters, dt)


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Generate cells with given hematocrit in pipe.')
    parser.add_argument('--params', type=str, default="parameters.pkl", help="Parameters file.")
    parser.add_argument('--abf-coords', type=str, required=True, help="The coordinates of the frozen particles of the abf.")
    parser.add_argument('--ranks', type=int, nargs=3, default=(1, 1, 1), help="Ranks in each dimension.")
    parser.add_argument('--Ht', type=float, default=0.20, help="Hematocrit, in [0,1].")
    parser.add_argument('--scale-ini', type=float, default=0.5, help="Initial size of the RBCs.")
    parser.add_argument('--seed', type=int, default=12345, help="Seed for initial placement of the RBCs.")
    args = parser.parse_args(argv)

    abf_coords = np.loadtxt(args.abf_coords)

    p = load_parameters(args.params)

    path = np.load(p.path_filename)
    generate_cells(p,
                   abf_coords=abf_coords.tolist(),
                   Ht=args.Ht,
                   scale_ini=args.scale_ini,
                   ranks=tuple(args.ranks),
                   seed=args.seed,
                   position_abf = list(path[0]))

if __name__ == '__main__':
    main(sys.argv[1:])
