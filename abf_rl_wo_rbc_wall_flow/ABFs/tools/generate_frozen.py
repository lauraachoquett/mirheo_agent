#!/usr/bin/env python3

import numpy as np

from mpi4py import MPI

import dpdprops

def create_from_mesh(num_density: float,
                     vertices: np.ndarray,
                     triangles: np.ndarray,
                     inertia: np.ndarray,
                     comm: MPI.Comm,
                     niter: int=1000):

    import mirheo as mir

    def recenter(coords, com):
        coords = [[r[0]-com[0], r[1]-com[1], r[2]-com[2]] for r in coords]
        return coords

    # arbitrary
    nu = 20
    sound_speed = 20
    rc = 1

    p = dpdprops.create_dpd_params_from_props(kinematic_viscosity=nu, sound_speed=20, nd=num_density, rc=rc)

    dt = p.get_max_dt() / 2

    # bounding box
    bb_hi = np.array(vertices).max(axis=0).tolist()
    bb_lo = np.array(vertices).min(axis=0).tolist()

    mesh_size = [hi-lo for lo, hi in zip(bb_lo, bb_hi)]

    fact = 2.0
    domain = (float (int (fact * mesh_size[0] + 1) ),
              float (int (fact * mesh_size[1] + 1) ),
              float (int (fact * mesh_size[2] + 1) ))

    ranks  = (1, 1, 1)

    u = mir.Mirheo(ranks, domain, debug_level=0, log_filename='', no_splash=True, comm_ptr=MPI._addressof(comm))

    dpd =  mir.Interactions.Pairwise('dpd', rc, kind="DPD", **p.to_interactions())
    vv = mir.Integrators.VelocityVerlet('vv')

    coords = [bb_lo, bb_hi]
    com_q  = [[0.5 * domain[0], 0.5 * domain[1], 0.5 * domain[2],   1., 0, 0, 0]]

    mesh = mir.ParticleVectors.Mesh(vertices, triangles)

    fakeOV = mir.ParticleVectors.RigidObjectVector('OV', mass=p.mass, inertia=tuple(inertia), object_size=len(coords), mesh=mesh)
    fakeIc = mir.InitialConditions.Rigid(com_q, coords)
    belongingChecker = mir.BelongingCheckers.Mesh("meshChecker")

    pvMesh = u.makeFrozenRigidParticles(belongingChecker, fakeOV, fakeIc, [dpd], vv, num_density, niter, dt)

    if pvMesh:
        frozenCoords = pvMesh.getCoordinates()
        frozenCoords = recenter(frozenCoords, com_q[0])
    else:
        frozenCoords = [[]]

    if u.isMasterTask():
        return frozenCoords
    else:
        return None

def make_positive(v):
    i = np.argmax(abs(v))
    if v[i] < 0:
        return -v
    return v

def shift_mesh(mesh):
    mesh.vertices -= mesh.center_mass
    return mesh

def rotate_mesh(mesh):
    return mesh
    # w, v = np.linalg.eig(mesh.moment_inertia)
    # v[0] = make_positive(v[0])
    # v[1] = make_positive(v[1])
    # v[2] = make_positive(v[2])
    # R = np.transpose(v)
    # mesh.vertices = np.transpose(np.matmul(R, mesh.vertices.T))
    # return mesh

def scale_mesh_to_length(mesh, length: float):
    extents = np.ptp(mesh.vertices, axis=0)
    assert len(extents) == 3
    L = np.max(extents)
    scale = length / L
    mesh.vertices = scale * mesh.vertices
    return mesh

def create_from_mesh_file(num_density: float,
                          mesh: str,
                          out_mesh: str,
                          length: float,
                          comm: MPI.Comm):
    import trimesh
    m = rotate_mesh(shift_mesh(trimesh.load(mesh)))
    if length is not None:
        m = scale_mesh_to_length(m, length)
    inertia = np.diag(m.moment_inertia)
    m.export(out_mesh)
    return create_from_mesh(num_density, m.vertices.tolist(), m.faces.tolist(), inertia, comm=comm)

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser("Freeze DPD particles inside a mesh, optionally rescaled.")
    parser.add_argument('mesh', type=str, help="Input mesh, will be rescaled.")
    parser.add_argument('--num-density', type=float, default=10,   help="Particle number density.")
    parser.add_argument('--length',      type=float, default=None, help="If set, the mesh will be rescaled to this length.")
    parser.add_argument('--out-coords',  type=str,   default="coords.txt",   help="File name (.txt) that will contain the coordinates of the frozen particles.")
    parser.add_argument('--out-mesh',    type=str,   default="ABF_mesh.ply", help="File name (.ply) that will contain the rescaled mesh.")
    args = parser.parse_args()

    rank = MPI.COMM_WORLD.rank
    key = rank // 2
    comm = MPI.COMM_WORLD.Split(color=rank, key=key)

    if key == 0:
        coords = create_from_mesh_file(args.num_density, args.mesh, args.out_mesh, args.length, comm=comm)

        if coords is not None:
            np.savetxt(args.out_coords, coords)

    MPI.COMM_WORLD.barrier()
