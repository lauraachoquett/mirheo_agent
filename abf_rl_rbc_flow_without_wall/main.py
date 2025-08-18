#!/usr/bin/env python

try:
    import cupy as cp
except Exception as e:
    print("cupy is not installed. This will cause an error if control is set.")
    cp = None
import torch

from utils import coordinate_in_path_ref_3D,coordinate_in_global_ref_3D
import time

import quaternion
import mirheo as mir
import numpy as np
import sys
from mpi4py import MPI
from utils import generate_simple_line,double_reflection_rmf,generate_helix,paraview_export,resample_path
from parameters import (Parameters,
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



def run_capillary_flow(p: 'Parameters',
                       *,
                       abf_coords: list,
                       comm: MPI.Comm,
                       ranks: tuple=(1,1,1),
                       restart_directory: str=None,
                       checkpoint_directory: str=None,
                        path_policy: str=None,
                        with_previous_action: bool = False,):

    rc = p.rc
    Lx = p.Lx
    Ly = p.Ly
    Lz = p.Lz
    RA = p.RA
    v = p.mean_vel
    max_contact_force = p.rbc_params.shear_modulus() * RA * 0.5
    print("With RBCs :",with_previous_action)
    dt = p.dpd_ii.get_max_dt()
    tend = 200 * 2 * np.pi / p.omega
    niters = int(tend/dt)

    domain = (Lx, Ly, Lz)

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

    vv = mir.Integrators.VelocityVerlet("vv")
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

    ## WALL :
    wall_0 = mir.Walls.MovingPlane(
        "wall_0",
        normal=[0.0, 0.0, -1.0],
        pointThrough=[0.0, 0.0, rc],
        velocity=[v, 0.0, 0.0],
    )
    wall_1 = mir.Walls.MovingPlane(
        "wall_1",
        normal=[0.0, 0.0, 1.0],
        pointThrough=[0.0, 0.0, Ly - rc ],
        velocity=[v, 0.0, 0.0],
    )

    u.registerWall(wall_0)
    u.registerWall(wall_1)

    pv_frozen_1 = u.makeFrozenWallParticles(
        pvName="wall_0",
        walls=[wall_0],
        interactions=[dpd_oo],
        integrator=vv,
        number_density=p.nd,
        dt=dt,
    )
    pv_frozen_2 = u.makeFrozenWallParticles(
        pvName="wall_1",
        walls=[wall_1],
        interactions=[dpd_oo],
        integrator=vv,
        number_density=p.nd,
        dt=dt,
    )

    u.setWall(wall_0, pv_outer)
    u.setWall(wall_1, pv_outer)

    # Set interactions between pvs

    u.registerInteraction(dpd_rbco)
    u.registerInteraction(dpd_rbci)

    u.setInteraction(dpd_oo, pv_outer, pv_outer)
    u.setInteraction(dpd_ii, pv_inner, pv_inner)
    u.setInteraction(dpd_io, pv_inner, pv_outer)

    u.setInteraction(dpd_rbco, pv_rbc, pv_outer)
    u.setInteraction(dpd_rbci, pv_rbc, pv_inner)

    u.setInteraction(dpd_oo, pv_outer, pv_frozen_1)
    u.setInteraction(dpd_oo, pv_outer, pv_frozen_1)
    u.setInteraction(dpd_io, pv_inner, pv_frozen_1)
    
    u.setInteraction(dpd_io, pv_inner, pv_frozen_2)
    u.setInteraction(dpd_rbco, pv_rbc, pv_frozen_2)
    u.setInteraction(dpd_rbco, pv_rbc, pv_frozen_2)

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
    u.registerPlugins(mir.Plugins.createDumpMesh("rbc_mesh_dump", pv_rbc, dump_every, 'ply/'))

    u.registerPlugins(mir.Plugins.createDumpMesh("abf_mesh_dump", pv_abf, dump_every, 'ply/'))
    u.registerPlugins(mir.Plugins.createStats('stats.csv', every=stats_every))
    u.registerPlugins(mir.Plugins.createDumpObjectStats("obj_stats", ov=pv_abf, dump_every=dump_every, filename="obj_stats/abf.csv"))

    
    
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
        dt_control = T_precession
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
                p_end = np.array([domain[0], domain[1]*1/2, domain[2]*1/2])
                p0 =  np.array([0, domain[1]*1/2, domain[2]*1/2])
                # LINE OR HELIX
                # path ,d = generate_simple_line(p0,p_end,100000)
                path = generate_helix(num_points=10000, radius=Ly*1/4, pitch=Lx, turns=1, clockwise=True, x_0=0, y_0=4/5*domain[1], z_0=1/2*domain[2])
                # Resample path if necessary : 
                dist = np.array([abs(path[i + 1] - path[i]) for i in range(len(path) - 1)])
                distance_between_points = 3.5e-2
                n = ceil(np.max(dist) / distance_between_points)
                if n > 1:
                    path, distances = resample_path(path, len(path) * n)
                    
                print("Distance between point path :",np.linalg.norm(path[1]-path[0]))
                T_rmf, N_rmf, B_rmf = double_reflection_rmf(path)
                
                if u.isMasterTask():
                    paraview_export(path, f'paraview_export/')
                
                
            start_time = time.time()
            if ((t-dt)%dt_control <= 1*dt and t > 1*dt_control) or init:
                ce = pv_abf.local.per_object['com_extents']
                num_abfs = ce.__cuda_array_interface__['shape'][0]
                has_abf   = (num_abfs > 0)
                if has_abf: ## Rank owner 
                    
                    ### Position of the ABF : 
                    
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
                    print(f"[Rank {comm.Get_rank()}] has ABF and will send info  : {info}")
                        
                    ### Save states
                    policy.history_state()
                        
            ## Share information : 
            if ((t-dt)%dt_control <= 1*dt and t > 1*dt_control) or init:
                abf_rank.share_information(info,compute_comm_bis)
                
                # print(f"[Rank {comm.Get_rank()}] received this information : previous_x : {abf_rank.previous_x} and previous_action : {abf_rank.previous_action}")
                # print('-----------------------------------------------------------')


            omega_t = omega * t
            ex = (1, 0, 0)
            ## Compute the direction with the last action compute by the rank owner
            B_dir = coordinate_in_global_ref_3D(np.zeros(3), abf_rank.previous_action,abf_rank.previous_t,abf_rank.previous_n,abf_rank.previous_b)

            q = get_quaternion_between_vectors(ex, B_dir)
            q = quaternion.from_float_array(q)
            B = (0.0,
                B_magn * np.cos(omega_t),
                B_magn * np.sin(omega_t))

            output = quaternion.rotate_vectors(q, B)
            
            if ((t-dt)%dt_control <= 1*dt and t > 1*dt_control) or init : ## Print ##
                print(f"[Rank {comm.Get_rank()}] has direction : {B_dir} and magnetic field : {output} (B : {B})", flush=True)
                print(f"[Rank {comm.Get_rank()}] Time : {(time.time() - start_time) * 1e3} ms")

                init = False
                
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
        print("time between update action: ",dt_control)

        sys.stdout.flush()

    if restart_directory is not None:
        u.restart(restart_directory)

    u.dumpWalls2XDMF([wall_0], h=(1,1,1), filename="h5/wall_0")
    u.dumpWalls2XDMF([wall_1], h=(1,1,1), filename="h5/wall_1")

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
    parser.add_argument('--with_previous_action', action='store_true', default=False, help="With previous action")
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
        with_previous_action = args.with_previous_action
                       )

if __name__ == '__main__':
    main(sys.argv[1:])
