#!/usr/bin/env python

import mirheo as mir
import numpy as np
import sys

from parameters import (SimulationParams,
                        load_parameters)
import objplacement

def gen_com_q_pipe(*,
                   abf_position: np.ndarray,
                   domain: tuple,
                   L: float,
                   R: float,
                   rbc_mesh,
                   abf_mesh,
                   scale_ini: float,
                   Ht: float,
                   seed: int=4242):
    """
    Create initial positions for rbcs in a pipe, scaled to a smaller mesh.
    """
    rbc_extents = np.ptp(rbc_mesh.vertices, axis=0) * scale_ini
    abf_extents = np.ptp(abf_mesh.vertices, axis=0)
    rbc_radius = np.max(rbc_extents) / 2

    base_geom = objplacement.Cylinder(L=domain,
                                      radius=R,
                                      center=np.array(domain)/2,
                                      axis=0)

    abf_positions = abf_position.reshape((-1,3))

    geom = objplacement.GeometryWithHoles(base=base_geom,
                                          hole_positions=abf_positions,
                                          hole_extents=abf_extents)

    com_q = objplacement.generate_ic(geometry=geom,
                                     obj_volume=rbc_mesh.volume,
                                     obj_extents=rbc_extents,
                                     obj_radius=rbc_radius,
                                     target_volume_fraction=Ht,
                                     seed=seed)

    return com_q


def generate_cells(p,
                   *,
                   abf_coords: list,
                   Ht: float,
                   abf_uncenter_factor: float,
                   ranks: tuple=(1,1,1),
                   seed: int=4242):

    rc = p.rc
    L  = p.L
    R  = p.R
    RA = p.RA
    kb = p.rbc_params.bending_modulus()

    assert abf_uncenter_factor <= 1
    assert abf_uncenter_factor >= 0

    max_contact_force = p.rbc_params.shear_modulus() * RA * 10
    drag = 2

    domain = (L, 2*R + 4*rc, 2*R + 4*rc)

    abf_radius = p.abf_radius
    abf_max_radial_pos = R - abf_radius * 1.1 # safety
    abf_position = np.array([domain[0]/5,1/2 * domain[1], 1/2* domain[2]])
    abf_position[1] -= abf_uncenter_factor * abf_max_radial_pos



    T = np.sqrt(p.m * RA**2 / kb)
    tend = 2 * T
    dt = 0.001

    niters = int(tend/dt)
    checkpoint_every = niters-5

    u = mir.Mirheo(ranks, domain, debug_level=0, log_filename='generate_ic',
                   checkpoint_folder="generate_ic/", checkpoint_every=checkpoint_every)


    sigma = RA/10
    eps = 3 * kb

    contact = mir.Interactions.Pairwise('contact', rc, kind='RepulsiveLJ', epsilon=eps, sigma=sigma, aware_mode='Object', max_force=max_contact_force)
    u.registerInteraction(contact)



    vv = mir.Integrators.VelocityVerlet('vv')
    u.registerIntegrator(vv)



    mesh_abf = mir.ParticleVectors.Mesh(p.mesh_abf.vertices.tolist(), p.mesh_abf.faces.tolist())
    I = np.diag(p.mesh_abf.moment_inertia * p.nd * p.m).tolist()
    pv_abf = mir.ParticleVectors.RigidObjectVector('abf', mass=p.m, inertia=I, object_size=len(abf_coords), mesh=mesh_abf)
    ic_abf = mir.InitialConditions.Rigid([abf_position.tolist() + [1, 0, 0, 0]], abf_coords)
    # We fix th abf on purpose, no integrator for it.
    #vv_abf = mir.Integrators.RigidVelocityVerlet("vv_abf")
    u.registerParticleVector(pv_abf, ic_abf)



    if u.isMasterTask():
        print(f"domain={domain}")
        print(f"L={L}, R={R}, RA={RA}")
        sys.stdout.flush()


    dump_every = int(niters//5)
    u.registerPlugins(mir.Plugins.createDumpMesh("abf_mesh_dump", pv_abf, dump_every, 'ply_generate_ic/'))
    u.registerPlugins(mir.Plugins.createStats('stats', every=dump_every))

    u.run(niters, dt)
    print("RBCs and MS are generated")


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Generate cells with given hematocrit in pipe.')
    parser.add_argument('parameters', type=str, help="Parameters of the simulation.")
    parser.add_argument('--abf-coords', type=str, required=True, help="The coordinates of the frozen particles of the abf.")
    parser.add_argument('--Ht', type=float, default=0.106, help="Hematocrit, in [0,1].")
    parser.add_argument('--abf-uncenter-factor', type=float, default=0, help="How uncenter the abf is placed, in [0,1]: 0 = centerd, 1 = at the border of the pipe")
    parser.add_argument('--seed', type=int, default=4242, help="Random seed for initial positions.")
    parser.add_argument('--ranks',          type=int, nargs=3, default=(1, 1, 1), help="Ranks in each dimension.")
    args = parser.parse_args(argv)

    abf_coords = np.loadtxt(args.abf_coords)

    p = load_parameters(args.parameters)

    generate_cells(p,
                ranks=tuple(args.ranks),
                   abf_coords=abf_coords.tolist(),
                   Ht=args.Ht,
                   abf_uncenter_factor=args.abf_uncenter_factor,
                   seed=args.seed)

if __name__ == '__main__':
    main(sys.argv[1:])
