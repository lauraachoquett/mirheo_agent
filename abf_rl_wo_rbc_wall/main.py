#!/usr/bin/env python

try:
    import cupy as cp
except Exception as e:
    print("cupy is not installed. This will cause an error if control is set.")
    cp = None
import torch

import quaternion
import mirheo as mir
import numpy as np
import sys
from mpi4py import MPI
from utils import generate_simple_line,double_reflection_rmf,generate_helix,paraview_export
from parameters import (SimulationParams,
                        load_parameters)

from load_NN import (load_policy)

def get_quaternion_between_vectors(u, v):
    k_cos_theta = np.dot(u, v)
    k = np.sqrt(np.dot(u, u) * np.dot(v, v))

    if k_cos_theta == -k:
        return np.array([0, 0, 1, 0])

    q = np.array([k_cos_theta + k, *np.cross(u, v)])
    return q / np.linalg.norm(q)



def run_capillary_flow(p: 'SimulationParams',
                       *,
                       abf_coords: list,
                       comm: MPI.Comm,
                       ranks: tuple=(1,1,1),
                       restart_directory: str=None,
                       checkpoint_directory: str=None,
                        path_policy: str=None,
                        with_rbc: bool = False,):

    rc = p.rc
    L  = p.L
    R  = p.R
    RA = p.RA
    f = p.force
    f=0
    max_contact_force = p.rbc_params.shear_modulus() * RA * 0.5
    print("With RBCs :",with_rbc)
    dt = p.dpd_ii.get_max_dt()
    tend = 200 * 2 * np.pi / p.omega
    niters = int(tend/dt)

    domain = (L, 2*R + 4*rc, 2*R + 4*rc)

    # uargs = {
    #     'nranks': ranks,
    #     'domain': domain,
    #     'debug_level': 0,
    #     'log_filename': 'log'
    # }
    
    uargs = {'nranks': ranks,
             'domain': domain,
             'debug_level': 3,
             'log_filename': 'log',
             'comm_ptr': MPI._addressof(comm),
             'max_obj_half_length': RA
    }
    print("Uargs :",uargs)
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

    # ABF

    mesh_abf = mir.ParticleVectors.Mesh(p.mesh_abf.vertices.tolist(), p.mesh_abf.faces.tolist())
    I = np.diag(p.mesh_abf.moment_inertia * p.nd * p.m).tolist()
    pv_abf = mir.ParticleVectors.RigidObjectVector('abf', mass=p.m, inertia=I, object_size=len(abf_coords), mesh=mesh_abf)
    ic_abf = mir.InitialConditions.Restart('generate_ic')

    u.registerParticleVector(pv_abf, ic_abf)

    vv_abf = mir.Integrators.RigidVelocityVerlet("vv_abf")
    u.registerIntegrator(vv_abf)


    # Solvent

    pv_outer = mir.ParticleVectors.ParticleVector('outer', mass=p.m)
    ic_outer = mir.InitialConditions.Uniform(number_density=p.nd)
    u.registerParticleVector(pv=pv_outer, ic=ic_outer)

    rbc_checker = mir.BelongingCheckers.Mesh('rbc_checker')
    u.registerObjectBelongingChecker(rbc_checker, pv_rbc)
    pv_inner = u.applyObjectBelongingChecker(rbc_checker, pv_outer, inside='inner', correct_every=10000)

    # Interactions

    dpd_oo = mir.Interactions.Pairwise('dpd_oo', rc=rc, kind="DPD", **p.dpd_oo.to_interactions())
    dpd_ii = mir.Interactions.Pairwise('dpd_ii', rc=rc, kind="DPD", **p.dpd_ii.to_interactions())
    dpd_io = mir.Interactions.Pairwise('dpd_io', rc=rc, kind="DPD", **p.dpd_io.to_interactions())

    u.registerInteraction(dpd_oo)
    u.registerInteraction(dpd_ii)
    u.registerInteraction(dpd_io)

    dpd_rbco = mir.Interactions.Pairwise('dpd_rbco', rc=rc, kind="DPD", **p.dpd_rbco.to_interactions())
    dpd_rbci = mir.Interactions.Pairwise('dpd_rbci', rc=rc, kind="DPD", **p.dpd_rbci.to_interactions())

    rbc_int = mir.Interactions.MembraneForces('int_rbc', **p.rbc_params.to_interactions_zero_visc(), stress_free=True)
    u.registerInteraction(rbc_int)

    sigma = RA/7
    eps = 5 * p.rbc_params.bending_modulus()
    contact = mir.Interactions.Pairwise('contact', rc, kind='RepulsiveLJ', epsilon=eps, sigma=sigma, aware_mode='Object', max_force=max_contact_force)
    u.registerInteraction(contact)

    # Walls

    wall = mir.Walls.Cylinder('pipe', center=(domain[1]/2, domain[2]/2),
                              radius=R, axis='x', inside=True)
    u.registerWall(wall)
    frozen = u.makeFrozenWallParticles("frozen", walls=[wall], interactions=[dpd_oo], integrator=vv, number_density=p.nd, nsteps=int(1.0/dt), dt=dt)
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

    u.setInteraction(dpd_oo, pv_outer, pv_abf)
    u.setInteraction(dpd_io, pv_inner, pv_abf)
    u.setInteraction(dpd_rbco, pv_rbc, pv_abf)

    u.setInteraction(contact, pv_rbc, pv_rbc)
    u.setInteraction(contact, pv_rbc, pv_abf)
    u.setInteraction(rbc_int, pv_rbc, pv_rbc)

    # Integrators

    dt_rbc_el = p.rbc_params.get_max_dt_elastic(mass=p.m)
    dt_rbc_visc = p.rbc_params.get_max_dt_visc(mass=p.m)
    substeps = 5 + int(dt/dt_rbc_el)

    vv_rbc = mir.Integrators.SubStepShardlowSweep("vv_rbc", substeps=substeps,
                                                  fastForces=rbc_int, **p.rbc_params.to_viscous(),
                                                  nsweeps=10)
    u.registerIntegrator(vv_rbc)

    u.setIntegrator(vv_rbc, pv_rbc)
    u.setIntegrator(vv, pv_outer)
    u.setIntegrator(vv, pv_inner)
    u.setIntegrator(vv_abf, pv_abf)

    # bouncers

    abf_bouncer = mir.Bouncers.Mesh("abf_bouncer", "bounce_maxwell", kBT=0.0)
    u.registerBouncer(abf_bouncer)

    rbc_bouncer = mir.Bouncers.Mesh("rbc_bouncer", "bounce_maxwell", kBT=0.0)
    u.registerBouncer(rbc_bouncer)

    u.setBouncer(abf_bouncer, pv_abf, pv_outer)
    u.setBouncer(rbc_bouncer, pv_rbc, pv_outer)
    u.setBouncer(rbc_bouncer, pv_rbc, pv_inner)

    belonging_checker = mir.BelongingCheckers.Mesh("inside_abf_checker")
    u.registerObjectBelongingChecker(belonging_checker, pv_abf)
    u.applyObjectBelongingChecker(belonging_checker, pv_outer, correct_every=10000, inside="none", outside="")


    # Plugins

    t_dump_every = 0.01 * 2 * np.pi / p.omega

    stats_every = int(t_dump_every/dt)
    dump_every = int(t_dump_every/dt)   
    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force_rbc", pv_rbc, wall, C=1, h=rc, max_force=max_contact_force))
    u.registerPlugins(mir.Plugins.createDumpMesh("rbc_mesh_dump", pv_rbc, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force_abf", pv_abf, wall, C=10*p.dpd_oo.a/rc, h=rc, max_force=10*p.dpd_oo.a))

    u.registerPlugins(mir.Plugins.createDumpMesh("abf_mesh_dump", pv_abf, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createStats('stats.csv', every=stats_every))
    u.registerPlugins(mir.Plugins.createDumpObjectStats("obj_stats", ov=pv_abf, dump_every=dump_every, filename="obj_stats/abf.csv"))

    
    
    if path_policy is not None:
        import cupy as cp
        import quaternion
        compute_comm = comm.Split(color=int(u.isComputeTask()), key=comm.Get_rank())
        ngpus = cp.cuda.runtime.getDeviceCount()

        mygpu = compute_comm.rank % ngpus
        torch.cuda.set_device(mygpu)
        cp.cuda.runtime.setDevice(mygpu)

        policy = load_policy(path_policy, mygpu, 5,False)
        B_magn = p.magn_B
        m_magn = p.magn_m
        omega = p.omega

        state = u.getState()
        T_precession = 2 * np.pi / omega
        dt_control = 2 * T_precession
        control_update_every = int(dt_control / dt) 
        # Variables to store initialization state and path data
        init = True
        path = None
        T_rmf = N_rmf = B_rmf = None
        B_dir = None  # Start as None to detect first computation
        
        def magnetic_field_3D(t):
            nonlocal init, path, T_rmf, N_rmf, B_rmf, B_dir,t_dump_every
            p_end = np.array([domain[0], domain[1]*1/2, domain[2]*1/2])
            new_action=False
            if init:
                ce = pv_abf.local.per_object['com_extents']
                com_extents_abf = cp.asarray(ce)
                pos = com_extents_abf[0, 0:3].get()
                pos = state.domain_info.local_to_global(pos)
                print("INIT : pos abf ",pos)
                r0 = [pos[0], pos[1], pos[2]]
                p0 = np.array(r0)
                # INIT PATH
                # path ,d = generate_simple_line(p0,p_end,100000)
                path = generate_helix(num_points=50000, radius=R*1/3, pitch=130, turns=1, clockwise=True, x_0=pos[0], y_0=pos[1], z_0=pos[2])
                print("Distance between point path :",np.linalg.norm(path[1]-path[0]))
                T_rmf, N_rmf, B_rmf = double_reflection_rmf(path)
                paraview_export(path, 'paraview_export/')
                B_dir = np.array(T_rmf[0])
                print("INIT : B ",B_dir)
                policy.init_state(p0, path, T_rmf, N_rmf, B_rmf)
                init = False  # Set init to False to avoid recomputation

            else : 
                if ((t-dt)%dt_control <= 1*dt and t > 1*dt_control):
                    ce = pv_abf.local.per_object['com_extents']
                    com_extents_abf = cp.asarray(ce)
                    pos = com_extents_abf[0, 0:3].get()
                    pos = state.domain_info.local_to_global(pos)
                    r = np.array([pos[0], pos[1], pos[2]])
                    new_action = True
                    print("New Action")
                    policy.compute_state(r, path, T_rmf, N_rmf, B_rmf, dt*control_update_every)
                B_dir = np.array(policy(new_action))
                if new_action :
                    print('Direction : ',B_dir)
                    print('-----------------------------------------------------------')
                    
            B_dir /= np.linalg.norm(B_dir)
                        
            omega_t = omega * t
            ex = (1, 0, 0)

            q = get_quaternion_between_vectors(ex, B_dir)
            q = quaternion.from_float_array(q)
            B = (0.0,
                 B_magn * np.cos(omega_t),
                 B_magn * np.sin(omega_t))
            
            if (t%(dump_every*dt) <= 1*dt):
                policy.history_state()
            output = quaternion.rotate_vectors(q, B)
            return output
        
        magnetic_field = magnetic_field_3D
         
    u.registerPlugins(mir.Plugins.createExternalMagneticTorque("magnetic_torque", pv_abf, m_magn, magnetic_field))

    if u.isMasterTask():
        print(f"tend = {tend}")
        print(f"dt = {dt}")
        print(f"substeps = {substeps}")
        print(f"omega = {omega}")
        print(f"dump every = {dump_every}")
        print("control_update_every : ", control_update_every)
        print("time between update : ",dt_control)

        sys.stdout.flush()

    if restart_directory is not None:
        u.restart(restart_directory)

    u.run(niters, dt)


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Run rbcs flowing in a capillary pipe.')
    parser.add_argument('--parameters',         type=str, default="parameters.pkl", help="Parameters file.")
    parser.add_argument('--abf-coords', type=str, required=True, help="The coordinates of the frozen particles of the abf.")
    parser.add_argument('--control-param', type=float, default=0, help="The parameter used to control the abf trajectory. 0 to disable control.")
    parser.add_argument('--restart-dir', type=str, default=None, help="The restart directory name (no restart if not set)")
    parser.add_argument('--checkpoint-dir', type=str, default=None, help="The checkpoint directory name (no checkpoint if not set)")
    parser.add_argument('--no-visc', action='store_true', default=False, help="Disable viscosity")
    parser.add_argument('--with_rbc', action='store_true', default=False, help="With RBC")
    parser.add_argument('--policy', type=str, default=None, help="If set, the json file that contains the NN (from korali)")    
    parser.add_argument('--ranks',          type=int, nargs=3, default=(1, 1, 1), help="Ranks in each dimension.")
    
    args = parser.parse_args(argv)

    p = load_parameters(args.parameters)
    abf_coords = np.loadtxt(args.abf_coords)
    comm = MPI.COMM_WORLD

    run_capillary_flow(
        p, ranks=tuple(args.ranks),
        comm=comm,
        abf_coords=abf_coords,
        restart_directory=args.restart_dir,
        checkpoint_directory=args.checkpoint_dir,
        path_policy=args.policy,
        with_rbc = args.with_rbc
                       )

if __name__ == '__main__':
    main(sys.argv[1:])
