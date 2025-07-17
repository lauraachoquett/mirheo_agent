#!/usr/bin/env python

import mirheo as mir
import numpy as np
import os
import trimesh
import sys

import dpdprops

def run(*,
        mesh: "trimesh.Trimesh",
        vel: list,
        omega: list,
        coords_file: str,
        Re: float=0.1,
        Ma: float=0.05,
        dry_run: bool=False):

    extents = np.ptp(mesh.vertices, axis=0)
    assert len(extents) == 3

    body_length = max(extents)
    L = 5 * body_length # TODO

    ranks = (1, 1, 1)
    domain = (L, L, L)

    nd = 10
    m = 1
    rc = 1

    U = max(vel)
    if U == 0:
        U = body_length/2 * max(omega)

    dpd_prms = dpdprops.create_dpd_params_from_Re_Ma(Re=Re, Ma=Ma, U=U, L=body_length,
                                                     nd=nd, rc=rc, mass=m)

    u = mir.Mirheo(ranks, domain, debug_level=3, log_filename='log')

    pv_solvent = mir.ParticleVectors.ParticleVector('solvent', mass=m)
    ic_solvent = mir.InitialConditions.Uniform(number_density=nd)
    u.registerParticleVector(pv=pv_solvent, ic=ic_solvent)

    coords = np.loadtxt(coords_file).tolist()
    I = np.diag(mesh.moment_inertia * nd * m).tolist()
    rigid_mesh = mir.ParticleVectors.Mesh(mesh.vertices.tolist(), mesh.faces.tolist())
    com_q = [[domain[0]/2, domain[1]/2, domain[2]/2, 1., 0, 0, 0]]

    pv_rigid = mir.ParticleVectors.RigidObjectVector('rigid', m, inertia=I, object_size=len(coords), mesh=rigid_mesh)
    ic_rigid = mir.InitialConditions.Rigid(com_q, coords)
    u.registerParticleVector(pv=pv_rigid, ic=ic_rigid)

    rigid_checker = mir.BelongingCheckers.Mesh('rigid_checker')
    u.registerObjectBelongingChecker(rigid_checker, pv_rigid)
    u.applyObjectBelongingChecker(rigid_checker, pv_solvent, inside='none', correct_every=0)

    vv = mir.Integrators.VelocityVerlet("vv")
    vv_rigid = mir.Integrators.RigidVelocityVerlet("vv_rigid")
    u.registerIntegrator(vv)
    u.registerIntegrator(vv_rigid)
    u.setIntegrator(vv, pv_solvent)
    u.setIntegrator(vv_rigid, pv_rigid)

    dpd = mir.Interactions.Pairwise('dpd', rc, kind="DPD", **dpd_prms.to_interactions())
    u.registerInteraction(dpd)

    u.setInteraction(dpd, pv_solvent, pv_solvent)
    u.setInteraction(dpd, pv_rigid, pv_solvent)

    bouncer = mir.Bouncers.Mesh("bouncer", "bounce_maxwell", kBT=dpd_prms.kBT)
    u.registerBouncer(bouncer)
    u.setBouncer(bouncer, pv_rigid, pv_solvent)

    dt = dpd_prms.get_max_dt() / 2
    T = body_length / U
    tend = 2 * T

    stats_every = int(0.1 * T / dt)

    path = f"pin_U_{vel[0]}_{vel[1]}_{vel[2]}_W_{omega[0]}_{omega[1]}_{omega[2]}"

    V = [mir.Plugins.PinObject.Unrestricted if v < 0 else v for v in vel]
    W = [mir.Plugins.PinObject.Unrestricted if w < 0 else w for w in omega]

    u.registerPlugins(mir.Plugins.createPinObject('pin_object', pv_rigid, stats_every, path, velocity=V, angular_velocity=W))

    path = f"stats_U_{vel[0]}_{vel[1]}_{vel[2]}_W_{omega[0]}_{omega[1]}_{omega[2]}"
    u.registerPlugins(mir.Plugins.createDumpObjectStats('object_stats', pv_rigid, stats_every, path))

    u.registerPlugins(mir.Plugins.createStats("Stats", "stats.csv", stats_every))

    # set the velocity to zero far from the rigid object

    low = [0,0,0]
    high = [domain[0], domain[1], domain[2]]
    UW = list(vel) + list(omega)
    assert len(UW) == 6

    i = np.argmax(UW)
    if i < 3:
        d = (i+1)%3
        high[d] = rc
        u.registerPlugins(mir.Plugins.createImposeProfile("imposeVel", pv=pv_solvent, low=low, high=high, velocity=(0,0,0), kBT=dpd_prms.kBT))

    if u.isMasterTask():
        print(f"extents = {extents}")
        print(f"domain = {domain}")
        print(f"dpd parameters = {dpd_prms}")
        print(f"tend = {tend}")
        print(f"dt = {dt}")
        sys.stdout.flush()

    if not dry_run:
        u.run(int(tend/dt), dt)

        if u.isMasterTask():
            fname = os.path.join(path, "params.txt")
            with open(fname, "w") as f:
                print(f"eta = {dpd_prms.dynamic_viscosity()}", file=f)
                print(f"l = {body_length}", file=f)



def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Translate / rotate a rigid object and collect the hydrodynamic forces and torques.')
    parser.add_argument('mesh', type=str, help="Object mesh.")
    parser.add_argument('coords', type=str, help="Frozen particle locations.")
    parser.add_argument('--dry-run', default=False, action='store_true')
    parser.add_argument('--U', nargs=3, type=float, required=True, help="Linear velocity")
    parser.add_argument('--W', nargs=3, type=float, required=True, help="Angular velocity")
    args = parser.parse_args(argv)

    vel = args.U
    omega = args.W

    # check that there is at least one positive component
    assert max([np.max(vel), np.max(omega)]) > 0

    mesh = trimesh.load(args.mesh)

    run(mesh=mesh,
        coords_file=args.coords,
        vel=vel,
        omega=omega,
        dry_run=args.dry_run)

if __name__ == '__main__':
    main(sys.argv[1:])
