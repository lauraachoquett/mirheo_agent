#!/usr/bin/env python3

import json
import numpy as np
from TD3 import TD3
from utils import coordinate_in_path_ref_3D,coordinate_in_global_ref_3D
import torch
import matplotlib.pyplot as plt
import os 
import pyvista as pv
import pickle 
from mpi4py import MPI

class PolicyEvaluator:
    def __init__(self, policy_file: str, device_id: str, n_lookahead: int,  previous_action :bool,control_update_every:int,rank:int):
        self.device = torch.device(f"cuda:{device_id}" if torch.cuda.is_available() else "cpu")
        self.rank=rank
        self.n_lookahead = n_lookahead

        state_dim = 3 * (1 + 1 + n_lookahead)
        self.previous_action = previous_action
        if self.previous_action:
            state_dim+=3
            
        action_dim = 3
        max_action = np.inf
        self.agent = TD3(state_dim, action_dim, max_action)
        self.agent.load(policy_file, device=self.device)
        # print("Agent loaded on :", self.device)
        
        self.length_cylinder = 50
        self.typical_length_rom = 1
        self.length_scale = self.length_cylinder / self.typical_length_rom
        
        self.velocity_scale = 200 / 1.0 ## TODO : Rescale ? 
        
        self.previous_x = None
        self.x = np.zeros(3)
        self.past_action = None
        self.t = np.zeros(3)
        self.n = np.zeros(3)
        self.b = np.zeros(3)
        self.p_cp = None
        
        self.d=0
        self.beta=0.25
        self.C=1
        
        self.state = np.zeros(3*(1+1+1*self.n_lookahead))
        self.state_list={}
        
        self.x_list=[]



        
    def compute_state(self, x, path, T, N, B, dt,past_action,previous_x,init):
        if init :
            past_action = np.array(T[0])
            previous_x = x
        self.x_list.append(x)
        self.x = x
        self.previous_x = previous_x
        self.past_action = past_action
        try : 
            d, id_cp = min_dist_closest_point(self.x, path)
        except Exception as e:
            return self.state
        
        self.d=d
        path_len = len(path)
        self.p_cp = path[id_cp]
        self.t = T[id_cp]
        self.n = N[id_cp]
        self.b = B[id_cp]

        s = coordinate_in_path_ref_3D(self.p_cp, self.x, self.t, self.n, self.b)
        s_scaled = s/self.length_scale
        
        result = [s_scaled.reshape(1, 3)]
        
        s_previous = coordinate_in_path_ref_3D(self.p_cp, self.previous_x, self.t, self.n, self.b)
        v_local_path = ((s - s_previous) / dt)
        # if len(self.x_list)>=2:
        #     print("Velocity before scaling without length scaling:", (self.x_list[-1]-self.x_list[-2])/dt)
        #     print("Velocity before scaling with length scaling:", v_local_path)
        
        v_local_path *= self.velocity_scale
        v_norm = np.linalg.norm(v_local_path) * self.velocity_scale
        if v_norm > 2 and v_norm != 0:
            v_local_path = v_local_path / v_norm * 2
            
        result.append(v_local_path .reshape(1, 3))

        if self.n_lookahead > 0:
            lookahead = []
            for i in range(1, self.n_lookahead + 1):
                idx = min(id_cp + i, path_len - 1)
                next_p = path[idx]
                next_p_local = coordinate_in_path_ref_3D(self.p_cp, next_p, self.t, self.n, self.b).reshape(1, 3)
                lookahead.append(next_p_local)
                
            result.append(np.concatenate(lookahead, axis=0))
            
        if self.previous_action:
            
            result.append(self.past_action.reshape(1, 3))
            
        state_true = np.concatenate(result, axis=0)

        self.state= state_true

        self.previous_x = self.x
        
    def __call__(self,t):
        self.state_list[t]=self.state
        action = self.agent.select_action(self.state)
        self.past_action=action
        self.paraview_pos(f'paraview_export_rank_{self.rank}/') 

            


    def history_state(self):
        file_name = 'save_states/'
        os.makedirs(file_name, exist_ok=True)
        with open(os.path.join(file_name, "states.pkl"), "wb") as f:
            pickle.dump(self.state_list, f)
            

    def paraview_pos(self, output_save_path):
        os.makedirs(output_save_path, exist_ok=True)


        # Save x_list as PolyData
        if self.x_list is not None and len(self.x_list) > 0:
            path_polydata = pv.PolyData(self.x_list)
            if len(self.x_list) > 1:
                lines = []
                for i in range(len(self.x_list) - 1):
                    lines.extend([2, i, i + 1])
                path_polydata.lines = np.array(lines)
                print(f"Nombre de segments créés: {len(self.x_list) - 1}")
            path_output = os.path.join(output_save_path, 'true_position_abf.vtp')
            path_polydata.save(path_output)
            print(f"Chemin sauvegardé: {path_output}")

        

    def reward(self,x_target,dt):
        d = self.d
        rew_t = -self.C * dt
        rew_target = -np.linalg.norm(self.x - x_target) + np.linalg.norm(
            self.previous_x - x_target
        )
        rew_d = -self.beta * d
        rew = rew_t + rew_d + rew_target
        return rew_t, rew_d, rew_target, rew
    
    def forward_info(self):
        return [self.x,self.past_action,self.t,self.n,self.b]
    
def load_policy(policy_file: str, device_id: str, n_lookahead: int,previous_action:bool, control_update_every : int,rank:int) -> callable:
    return PolicyEvaluator(policy_file, device_id, n_lookahead,previous_action,control_update_every,rank)

def min_dist_closest_point(x, path):
    try:
        if x is None:
            raise ValueError(f"Input x is None : {x}")
        path = np.asarray(path)
        x = np.asarray(x)
        dists = np.linalg.norm(path - x, axis=-1)
        idx = np.argmin(dists)
        return dists[idx], idx
    except Exception as e:
        print(f"Error in min_dist_closest_point: {e}")
        raise



class ABFRankTracker:
    def __init__(self):
        self.previous_x = np.ones(3)*-np.inf
        self.previous_action = np.ones(3)*-np.inf
        self.previous_t = np.ones(3)*-np.inf
        self.previous_n = np.ones(3)*-np.inf
        self.previous_b = np.ones(3)*-np.inf
    
    def share_information(self, info, compute_comm_bis):
        if compute_comm_bis == MPI.COMM_NULL:
            return False  

        x_in = np.ascontiguousarray(info[0], dtype=np.float64)
        a_in = np.ascontiguousarray(info[1], dtype=np.float64)
        t_in = np.ascontiguousarray(info[2], dtype=np.float64)
        n_in = np.ascontiguousarray(info[3], dtype=np.float64)
        b_in = np.ascontiguousarray(info[4], dtype=np.float64)

        x_out = np.empty_like(x_in)
        a_out = np.empty_like(a_in)
        t_out = np.empty_like(t_in)
        n_out = np.empty_like(n_in)
        b_out = np.empty_like(b_in)

        compute_comm_bis.Allreduce(x_in, x_out, op=MPI.MAX)
        compute_comm_bis.Allreduce(a_in, a_out, op=MPI.MAX)
        compute_comm_bis.Allreduce(t_in, t_out, op=MPI.MAX)
        compute_comm_bis.Allreduce(n_in, n_out, op=MPI.MAX)
        compute_comm_bis.Allreduce(b_in, b_out, op=MPI.MAX)
        
        self.previous_x      = x_out
        self.previous_action = a_out
        self.previous_t = t_out
        self.previous_n = n_out
        self.previous_b = b_out
        
