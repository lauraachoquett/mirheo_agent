#!/usr/bin/env python
try:
    import cupy as cp
except Exception as e:
    print("cupy is not installed. This will cause an error if control is set.")
    cp = None
import torch

from utils import coordinate_in_path_ref_3D,coordinate_in_global_ref_3D
import time
import os, sys

import quaternion
import mirheo as mir
import numpy as np
import sys
from mpi4py import MPI
from utils import generate_simple_line,double_reflection_rmf,generate_helix,paraview_export,resample_path
from parameters import (Parameters,
                        ContactParams,
                        load_parameters)

from load_NN import (load_policy,ABFRankTracker)
from math import floor,ceil


def get_quaternion_between_vectors(u, v):
    k_cos_theta = np.dot(u, v)
    k = np.sqrt(np.dot(u, u) * np.dot(v, v))

    if k_cos_theta == -k:
        return np.array([0, 0, 1, 0])

    q = np.array([k_cos_theta + k, *np.cross(u, v)])
    return q / np.linalg.norm(q)

def run(p: Parameters,
        *,
        comm: MPI.Comm,
        abf_coords: list,
        no_visc: bool,
        ranks: tuple=(1,1,1),
        restart_directory: str=None,
        checkpoint_directory: str=None,
        path_policy: str=None
        ):
    """
    Arguments:
        p: the simulation parameters (see parameters.py)
        comm: MPI communicator
        abf_coords: list of corrdinates of frozen particles, in the abf frame of reference
        ranks: number of ranks per dimension
        restart_directory: if set, specify from which directory to restart.
        checkpoint_directory: if set, specify to which directory to dump checkpoints.
        path_policy: if None, no control; else, path to a json file containing a NN
    """

    rc     = p.rc
    domain = p.domain
    RA     = p.RA

    dt_fluid = min([dpd.get_max_dt() for dpd in [p.dpd_ii, p.dpd_oo]])

    tend = 1000 * RA / p.U

    uargs = {'nranks': ranks,
             'domain': domain,
             'debug_level': 3,
             'log_filename': 'log',
             'comm_ptr': MPI._addressof(comm),
             'max_obj_half_length': RA
    }

    if checkpoint_directory is not None:
        save_checkpoint_every = 10000
        uargs['checkpoint_every']  = save_checkpoint_every
        uargs['checkpoint_folder'] = checkpoint_directory

    u = mir.Mirheo(**uargs)

    vv = mir.Integrators.VelocityVerlet('vv')
    u.registerIntegrator(vv)

    # Membranes

    mesh_rbc = mir.ParticleVectors.MembraneMesh(p.mesh_ini.vertices.tolist(),
                                                p.mesh_ref.vertices.tolist(),
                                                p.mesh_ini.faces.tolist())
    pv_rbc = mir.ParticleVectors.MembraneVector('rbc', mass=p.m, mesh=mesh_rbc)
    ic_dir = os.path.abspath('generate_ic')
    
    if not os.path.isdir(ic_dir):
        print(f"[rank {comm.rank}] IC folder missing: {ic_dir}", file=sys.stderr)
    comm.Barrier()
    
    u.registerParticleVector(pv_rbc, mir.InitialConditions.Restart('generate_ic'))

    # abf

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
    pv_inner = u.applyObjectBelongingChecker(rbc_checker, pv_outer,
                                             inside='inner', correct_every=10000)

    # Interactions

    dpd_oo = mir.Interactions.Pairwise('dpd_oo', rc=rc, kind="DPD", **p.dpd_oo.to_interactions())
    dpd_ii = mir.Interactions.Pairwise('dpd_ii', rc=rc, kind="DPD", **p.dpd_ii.to_interactions())
    dpd_io = mir.Interactions.Pairwise('dpd_io', rc=rc, kind="DPD", **p.dpd_io.to_interactions())

    u.registerInteraction(dpd_oo)
    u.registerInteraction(dpd_ii)
    u.registerInteraction(dpd_io)

    dpd_rbco = mir.Interactions.Pairwise('dpd_rbco', rc=rc, kind="DPD",
                                         **p.dpd_rbco.to_interactions())
    dpd_rbci = mir.Interactions.Pairwise('dpd_rbci', rc=rc, kind="DPD",
                                         **p.dpd_rbci.to_interactions())

    # zero viscosity here, we will either use the Shardlow integrator or not viscosity at all.
    rbc_int = mir.Interactions.MembraneForces('int_rbc',
                                              **p.rbc_params.to_interactions_zero_visc(),
                                              stress_free=True)
    u.registerInteraction(rbc_int)

    cp = p.contact_rbcs
    sigma             = cp.sigma
    eps               = cp.eps
    max_contact_force = cp.max_contact_force

    contact = mir.Interactions.Pairwise('contact', rc, kind='RepulsiveLJ',
                                        epsilon=eps, sigma=sigma, aware_mode='Object',
                                        max_force=max_contact_force)
    u.registerInteraction(contact)

    # Walls

    wall = mir.Walls.SDF('walls', p.sdf_filename)
    u.registerWall(wall)
    frozen = u.makeFrozenWallParticles("frozen", walls=[wall], interactions=[dpd_oo],
                                       integrator=vv, number_density=p.nd,
                                       nsteps=1000, dt=p.dpd_oo.get_max_dt())
    u.setWall(wall, pv_outer)

    # Set interactions between pvs

    u.registerInteraction(dpd_rbco)
    u.registerInteraction(dpd_rbci)

    u.setInteraction(dpd_oo, pv_outer, pv_outer)
    u.setInteraction(dpd_ii, pv_inner, pv_inner)
    u.setInteraction(dpd_io, pv_inner, pv_outer)

    u.setInteraction(dpd_rbco, pv_rbc, pv_outer)
    u.setInteraction(dpd_rbci, pv_rbc, pv_inner)

    u.setInteraction(dpd_oo, pv_outer, pv_abf)
    u.setInteraction(dpd_io, pv_inner, pv_abf)
    u.setInteraction(dpd_rbco, pv_rbc, pv_abf)

    u.setInteraction(dpd_oo, pv_outer, frozen)
    u.setInteraction(dpd_io, pv_inner, frozen)
    u.setInteraction(dpd_rbco, pv_rbc, frozen)
    u.setInteraction(dpd_oo,   pv_abf, frozen)

    u.setInteraction(contact, pv_rbc, pv_rbc)
    u.setInteraction(contact, pv_rbc, pv_abf)

    # Bouncers

    rbc_bouncer = mir.Bouncers.Mesh("rbc_bouncer", "bounce_maxwell", kBT=p.kBT)
    u.registerBouncer(rbc_bouncer)

    abf_bouncer = mir.Bouncers.Mesh("abf_bouncer", "bounce_maxwell", kBT=p.kBT)
    u.registerBouncer(abf_bouncer)

    u.setBouncer(rbc_bouncer, pv_rbc, pv_inner)
    u.setBouncer(rbc_bouncer, pv_rbc, pv_outer)
    u.setBouncer(abf_bouncer, pv_abf, pv_outer)

    # Integrators

    u.setIntegrator(vv, pv_outer)
    u.setIntegrator(vv, pv_inner)

    dt_rbc_el = p.rbc_params.get_max_dt_elastic(mass=p.m)
    dt_rbc_visc = p.rbc_params.get_max_dt_visc(mass=p.m)
    substeps = 5 + int(dt_fluid / dt_rbc_el)
    dt = dt_fluid / 5

    if no_visc:
        vv_rbc = mir.Integrators.SubStep("vv_rbc", substeps=substeps,
                                         fastForces=[rbc_int])
    else:
        nsweeps = 10
        vv_rbc = mir.Integrators.SubStepShardlowSweep("vv_rbc", substeps=substeps,
                                                      fastForces=rbc_int,
                                                      **p.rbc_params.to_viscous(),
                                                      nsweeps=nsweeps)
    u.registerIntegrator(vv_rbc)
    u.setIntegrator(vv_rbc, pv_rbc)

    u.setIntegrator(vv_abf, pv_abf)

    # Plugins

    t_dump_every = 0.5 * RA / p.U
    dump_every = int(t_dump_every/dt)
    stats_every = dump_every

    h = p.wall_repulsion_length
    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force_rbc", pv_rbc, wall,
                                                      C=max_contact_force/h, h=h,
                                                      max_force=max_contact_force))
    u.registerPlugins(mir.Plugins.createWallRepulsion("wall_force_abf", pv_abf, wall,
                                                      C=max_contact_force/h, h=h,
                                                      max_force=max_contact_force))


    f = p.forces_filename
    h = (rc, rc, rc)
    u.registerPlugins(mir.Plugins.createAddForceField("body_force_solvent", pv_outer, f, h))
    u.registerPlugins(mir.Plugins.createAddForceField("body_force_hemoglobine", pv_inner, f, h))
    u.registerPlugins(mir.Plugins.createAddForceField("body_force_membranes", pv_rbc, f, h))

    u.registerPlugins(mir.Plugins.createDumpMesh("rbc_mesh_dump", pv_rbc, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createDumpMesh("abf_mesh_dump", pv_abf, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createStats("stats", every=stats_every,
                                              filename='stats.csv',
                                              pvs=[pv_inner, pv_outer, pv_rbc]))
    u.registerPlugins(mir.Plugins.createDumpObjectStats("abf_stats", ov=pv_abf,
                                                        dump_every=stats_every,
                                                        filename="obj_stats/abf.csv"))


    compute_comm = None
    mygpu = None
    policy = None
    
    B_magn = p.magn_B
    m_magn = p.magn_m
    omega = p.omega
    
    print(f"BEFORE THE LOOP", flush=True)
    rank = comm.Get_rank()
    print(f"[Rank {rank}] isComputeTask() = {u.isComputeTask()}", flush=True)
    
    if path_policy is not None:
        print("enter the loop")
        compute_comm = comm.Split(color=int(u.isComputeTask()), key=rank)
        ngpus = cp.cuda.runtime.getDeviceCount()
        mygpu = compute_comm.rank % ngpus
        


        for i in range(ngpus):
            print(f"[Rank {rank}] PyTorch device {i}: {torch.cuda.get_device_name(i)}")
        
        print(f"[Rank {rank}] compute_comm.rank = {compute_comm.rank}, mygpu = {mygpu}")
        cp.cuda.runtime.setDevice(mygpu)
        print(f"[Rank {rank}] Assigned GPU {mygpu}")

        state = u.getState()
        T_precession = 2 * np.pi / omega
        dt_control = 1/4*T_precession
        control_update_every = int(dt_control / dt) 

        
        
        init = True
        path = None
        T_rmf = N_rmf = B_rmf = None
        B_dir = np.zeros(3)
        
        is_compute = u.isComputeTask()

        compute_comm_bis = comm.Split(color=1 if is_compute else MPI.UNDEFINED,
                                key=rank)
        if is_compute:
            print(f"[Rank {rank}] set torch device on GPU {mygpu}...")
            torch.cuda.set_device(mygpu)
            print(f"[Rank {rank}] Loading policy on GPU {mygpu}...")

            policy = load_policy(policy_file = path_policy, device_id = mygpu, n_lookahead = 5,previous_action =  with_previous_action, control_update_every = control_update_every,rank=rank)
            print(f"[Rank {rank}] Policy loaded successfully.")
            torch.cuda.synchronize()

            print(f"[Rank {rank}] compute_comm_bis.rank = {compute_comm_bis.rank}")
            
            abf_rank = ABFRankTracker()
            print(f"[Rank {comm.Get_rank()}] class abf_rank init : OK")
            
        def magnetic_field_3D(t):
            nonlocal init, path, T_rmf, N_rmf, B_rmf, B_dir,t_dump_every,abf_rank
            info = np.ones((5,3))*-np.inf #previous_x, previous_action,previous_t,previous_n,previous_b

            if init:
                # INIT PATH
                path = np.load(p.path_filename)
                
                # Resample path if necessary : 
                dist = np.array([abs(path[i + 1] - path[i]) for i in range(len(path) - 1)])
                distance_between_points = 3.5e-2
                n = ceil(np.max(dist) / distance_between_points)
                if n > 1:
                    path, distances = resample_path(path, len(path) * n)
                    
                print("Distance between point path :",np.linalg.norm(path[1]-path[0]))
                
                ## Generate local frame : (RMF algorithm)
                T_rmf, N_rmf, B_rmf = double_reflection_rmf(path)
                
                if u.isMasterTask():
                    paraview_export(path, f'paraview_export/')
                
            if has_abf and restart_directory is not None :
                init != abf_rank.load(restart_directory)
                
            start_time = time.time()
            if ((t-dt)%dt_control <= 1*dt and t > 1*dt_control) or init:
                ce = pv_abf.local.per_object['com_extents']
                num_abfs = ce.__cuda_array_interface__['shape'][0]
                has_abf   = (num_abfs > 0)
                if has_abf: ## Rank owner 
                    
                    ### Position of the abf : 
                    com_extents_abf = cp.asarray(ce)
                    pos = com_extents_abf[0, 0:3].get()
                    pos = state.domain_info.local_to_global(pos)
                    r = np.array([pos[0], pos[1], pos[2]]) 
                    
                    ### Compute state :
                    policy.compute_state(r, path, T_rmf, N_rmf, B_rmf, dt*control_update_every,abf_rank.previous_action,abf_rank.previous_x,init)
                    
                    ### New action
                    policy(t)
                
                    ### previous_x, previous_action,previous_t,previous_n,previous_b put in info
                    info = policy.forward_info()
                    print(f"[Rank {comm.Get_rank()}] has abf and will send info  : {info}")
                        
                    ### Save states
                    policy.history_state()
                        
            ## Share information : 
            abf_rank.share_information(info,compute_comm_bis)
                
            # print(f"[Rank {comm.Get_rank()}] received this information : previous_x : {abf_rank.previous_x} and previous_action : {abf_rank.previous_action}")
            # print('-----------------------------------------------------------')


            omega_t = omega * t
            ex = (1, 0, 0)
            
            ## Compute the direction with the last action computed by the rank owner
            B_dir = coordinate_in_global_ref_3D(np.zeros(3), abf_rank.previous_action,abf_rank.previous_t,abf_rank.previous_n,abf_rank.previous_b)

            q = get_quaternion_between_vectors(ex, B_dir)
            q = quaternion.from_float_array(q)
            B = (0.0,
                B_magn * np.cos(omega_t),
                B_magn * np.sin(omega_t))

            output = quaternion.rotate_vectors(q, B)
            
            ## Print ## : 
            if ((t-dt)%dt_control <= 1*dt and t > 1*dt_control) or init : 
                print(f"[Rank {comm.Get_rank()}] has direction : {B_dir} and magnetic field : {output} (B : {B})", flush=True)
                print(f"[Rank {comm.Get_rank()}] Time : {(time.time() - start_time) * 1e3} ms")

                init = False
            
            if has_abf and (t-dt)%(save_checkpoint_every*dt)<=1*dt:
                abf_rank.save(checkpoint_directory)
                
            return output
        
    
    magnetic_field = magnetic_field_3D
         
    u.registerPlugins(mir.Plugins.createExternalMagneticTorque("magnetic_torque", pv_abf, m_magn, magnetic_field))


    if u.isMasterTask():
        print(f"tend = {tend}")
        print(f"dt = {dt}")
        print(f"substeps = {substeps}")
        if not no_visc:
            print(f"nsweeps = {nsweeps}")
        sys.stdout.flush()

    if restart_directory is not None:
        u.restart(restart_directory)
    else:
        u.dumpWalls2XDMF([wall], h=(1,1,1), filename='h5/wall')

    u.run(int(tend/dt), dt)


def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--params',         type=str, default="parameters.pkl", help="Parameters file.")
    parser.add_argument('--abf-coords', type=str, required=True, help="The coordinates of the frozen particles of the abf.")
    parser.add_argument('--ranks',          type=int, nargs=3, default=(1, 1, 1), help="Ranks in each dimension.")
    parser.add_argument('--restart-dir',    type=str, default=None, help="The restart directory name (no restart if not set)")
    parser.add_argument('--checkpoint-dir', type=str, default=None, help="The checkpoint directory name (no checkpoint if not set)")
    parser.add_argument('--no-visc', action='store_true', default=False, help="Disable viscosity")
    parser.add_argument('--policy', type=str, default=None, help="If set, the json file that contains the NN (from korali)")
    args = parser.parse_args(argv)

    p = load_parameters(args.params)
    abf_coords = np.loadtxt(args.abf_coords)

    comm = MPI.COMM_WORLD

    run(p, ranks=tuple(args.ranks),
        comm=comm,
        abf_coords=abf_coords,
        no_visc=args.no_visc,
        restart_directory=args.restart_dir,
        checkpoint_directory=args.checkpoint_dir,
        path_policy=args.policy,
        )


if __name__ == '__main__':
    main(sys.argv[1:])
