#!/usr/bin/env python

import mirheo as mir

import numpy as np
import pandas as pd
import pint
import sys
import trimesh

from dpdprops import (
    KantorWLCRBCDefaultParams,
    JuelicherLimRBCDefaultParams,
    load_stress_free_mesh,
    load_equilibrium_mesh,
)

from tools import compute_diameters, compute_micro_beads_forces


def rescale_by_area(mesh, A):
    mesh.vertices *= np.sqrt(A / mesh.area)
    return mesh


def stretch_and_compute_diameters(
    *,
    ureg,
    comm_address,
    mesh_ini,
    mesh_ref,
    force_: float,
    params,
    RA: float = 1,
    dump: bool = False,
):
    """
    Parameters:
        ureg: the unit registry
        comm_address: the address of the MPI communicator
        mesh_ini: Initial mesh
        mesh_ref: stress free mesh
        force_: The applied stretching force, in [pN]
        params: One of the model parameters.
        RA: equivalent radius of the membrane in simulation units (sets the length scale)
        dump: if True, dump ply files
    Returns:
        D0, D1: The diameters of the RBC, in [um]
    """

    @ureg.wraps(None, ureg.dimensionless)
    def to_sim(a):
        return a

    force_ *= ureg.pN

    A0 = 4 * np.pi * RA**2
    mesh_ini = rescale_by_area(mesh_ini, A0)
    mesh_ref = rescale_by_area(mesh_ref, A0)

    ranks = (1, 1, 1)

    safety = 3
    mass = 1

    domain = np.ptp(mesh_ini.vertices, axis=0) * safety
    domain = tuple(np.array(domain, dtype=int))

    u = mir.Mirheo(
        ranks,
        domain,
        debug_level=0,
        log_filename="log",
        no_splash=True,
        comm_ptr=comm_address,
    )

    mesh_rbc = mir.ParticleVectors.MembraneMesh(
        mesh_ini.vertices.tolist(), mesh_ref.vertices.tolist(), mesh_ini.faces.tolist()
    )
    pv_rbc = mir.ParticleVectors.MembraneVector("rbc", mass=mass, mesh=mesh_rbc)
    ic_rbc = mir.InitialConditions.Membrane(
        [[domain[0] * 0.5, domain[1] * 0.5, domain[2] * 0.5, 1.0, 0.0, 0.0, 0.0]]
    )
    u.registerParticleVector(pv_rbc, ic_rbc)

    force_scale = 0.25 * ureg.pN  # arbitrary
    length_scale = np.sqrt(params.A0 / A0)
    time_scale = 1 * ureg.s  # arbitrary
    mass_scale = (force_scale / length_scale * time_scale**2).to(ureg.kg)

    force = to_sim(force_ / force_scale)

    rbc_params = params.get_params(
        length_scale=length_scale,
        time_scale=time_scale,
        mass_scale=mass_scale,
        mesh=mesh_ini,
    )

    rbc_params.gamma = 200.0  # for killing oscillations

    tend = 50.0
    dt0 = rbc_params.get_max_dt(mass=mass)
    substeps = 500
    dt = substeps * dt0

    int_rbc = mir.Interactions.MembraneForces(
        "int_rbc", **rbc_params.to_interactions(), stress_free=True
    )
    u.registerInteraction(int_rbc)
    vv = mir.Integrators.SubStep("vv", substeps, [int_rbc])
    u.registerIntegrator(vv)
    u.setIntegrator(vv, pv_rbc)
    # u.setInteraction(int_rbc, pv_rbc, pv_rbc)

    dc_ = 2 * ureg.um
    dc = to_sim(dc_ / length_scale)

    forces = compute_micro_beads_forces(
        mesh_ini, contact_diameter=dc, bead_force=force
    ).tolist()

    u.registerPlugins(
        mir.Plugins.createMembraneExtraForce("stretchForce", pv_rbc, forces)
    )

    if dump:
        ply_path = f"ply_f_{to_sim(force_/ureg.pN)}"
        u.registerPlugins(
            mir.Plugins.createDumpMesh(
                "mesh_dump", pv_rbc, int(tend / 50 / dt), ply_path
            )
        )

    u.run(int(tend / dt), dt)

    is_master = u.isMasterTask()
    Dshort, Dlong = 0, 0

    if is_master:
        rbc_pos = pv_rbc.getCoordinates()
        Dlong, Dshort = compute_diameters(rbc_pos)

        # make it dimensional again
        Dlong = to_sim(Dlong * length_scale / ureg.um)
        Dshort = to_sim(Dshort * length_scale / ureg.um)
    del u

    return Dshort, Dlong, is_master


def main(argv):
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("mesh_ini", type=str, help="The initial mesh.")
    parser.add_argument("mesh_ref", type=str, help="The mesh of stress free state.")
    parser.add_argument(
        "exp_file", type=str, help="The csv file containing the experimental values."
    )
    parser.add_argument("--out", type=str, default="results.csv")
    parser.add_argument(
        "--dump",
        action="store_true",
        default=False,
        help="Will dump ply files if set to True.",
    )
    parser.add_argument(
        "--model", type=str, choices=["KantorWLC", "JuelicherLim"], default="KantorWLC"
    )
    args = parser.parse_args(argv)

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    comm_address = MPI._addressof(comm)

    ureg = pint.UnitRegistry()

    # see https://github.com/mikedh/trimesh/issues/338
    mesh_ini = trimesh.load_mesh(args.mesh_ini, process=False)
    mesh_ref = trimesh.load_mesh(args.mesh_ref, process=False)

    df = pd.read_csv(args.exp_file, sep=" ")

    D0s = list()
    D1s = list()

    if args.model == "KantorWLC":
        params = KantorWLCRBCDefaultParams(
            ureg, mu=4.058 * ureg.uN / ureg.m, kappab=5.906e-19 * ureg.J, x0=0.416
        )
    else:
        params = JuelicherLimRBCDefaultParams(ureg)
        mesh_ini = load_equilibrium_mesh()
        mesh_ref = load_stress_free_mesh()

    for force in df["force"]:
        D0, D1, is_master = stretch_and_compute_diameters(
            ureg=ureg,
            comm_address=comm_address,
            mesh_ini=mesh_ini,
            mesh_ref=mesh_ref,
            force_=force,
            params=params,
            dump=args.dump,
        )
        if is_master:
            print(force, D0, D1)
            sys.stdout.flush()

        D0s.append(D0)
        D1s.append(D1)

    df["D0sim"] = D0s
    df["D1sim"] = D1s

    if is_master:
        df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main(sys.argv[1:])
