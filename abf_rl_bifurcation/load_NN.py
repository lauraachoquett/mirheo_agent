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

class PolicyEvaluator:
    def __init__(self, policy_file: str, device_id: str, n_lookahead: int,  previous_action :bool):
        self.device = torch.device(f"cuda:{device_id}")
        self.n_lookahead = n_lookahead

        state_dim = 3 * (1 + 1 + n_lookahead)
        if previous_action:
            state_dim+=3
            
        action_dim = 3
        max_action = np.inf
        self.agent = TD3(state_dim, action_dim, max_action)
        self.agent.load(policy_file, device=self.device)

        self.radius_biggest_ca_p = 0.269
        self.radius_biggest_ca_grid = 6.43
        self.radius_biggest_phys = 20  # micrometer

        self.length_cylinder = 35
        self.typical_length_rom = 1
        self.length_scale = self.length_cylinder / self.typical_length_rom
        
        # self.time_rom = 3.333e-02 / 5
        # self.time_sim = 0.00163
        # self.time_scale = self.time_sim / self.time_rom
        
        self.velocity_scale =  360 / 1.0 ## TODO : Rescale ? 
        self.previous_x = None
        self.x = np.zeros(3)
        self.past_action = None
        self.t = np.zeros(3)
        self.n = np.zeros(3)
        self.b = np.zeros(3)
        
        self.d=0
        self.beta=0.25
        self.C=1
        
        self.state = np.zeros(3*(1+1+1*self.n_lookahead))
        self.state_list=[]
        
        self.pos_wo_window = []
        self.pos_wo_list=[]
        self.pos_wo_local_list = []
        
        self.x_list=[]
        self.x_list_local=[]

    def add_to_list_pos(self, pos):        
        self.pos_wo_window.append(coordinate_in_path_ref_3D(self.p_cp, pos, self.t, self.n, self.b))

    def reset_pos_wo(self):
        self.pos_wo_window = []

    def mean_pos_wo(self,pos):
        self.x_list.append(pos)
        x_local = coordinate_in_path_ref_3D(self.p_cp ,pos,self.t, self.n, self.b)
        self.x_list_local.append(x_local)
        
        pos_array = np.array(self.pos_wo_window)
        if len(pos_array) == 0:
            return None
        mean_pos_wo_local = np.mean(pos_array, axis=0)
        self.pos_wo_local_list.append(mean_pos_wo_local)
        
        return coordinate_in_global_ref_3D(self.p_cp,mean_pos_wo_local,self.t, self.n, self.b)
    
    def init_state(self, x, path, T, N, B):
        self.x_list.append(x)
        self.pos_wo_list.append(x)
        self.x = x
        self.previous_x = self.x
        self.past_action = np.array(T[0])
        path = path
        d, id_cp = min_dist_closest_point(self.x, path)
        self.d=d
        self.p_cp = path[id_cp]
        self.t = T[id_cp]
        self.n = N[id_cp]
        self.b = B[id_cp]
        
    def compute_state(self, x, path, T, N, B, dt):
        self.pos_wo_list.append(x)
        self.x = x
        
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
        v_local_path = ((s - s_previous) / dt)*self.velocity_scale
        # if len(self.x_list)>=2:
        #     print("Velocity before scaling without mean:", (self.x_list[-1]-self.x_list[-2])/dt)
        #     print("Velocity before scaling with mean:", v_local_path)
            
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
        # result.append(self.past_action.reshape(1, 3))
        state_true = np.concatenate(result, axis=0)
        ideal_state = np.array([s_scaled,
            [ 1.00000000e-01, -8.46208763e-08,  8.03898274e-07],
            [ 5.00005000e-05,  0.00000000e+00,  0.00000000e+00],
            [ 1.00001000e-04,  0.00000000e+00,  0.00000000e+00],
            [ 1.50001500e-04,  0.00000000e+00,  0.00000000e+00],
            [ 2.00002000e-04,  0.00000000e+00,  0.00000000e+00],
            [ 2.50002500e-04,  0.00000000e+00,  0.00000000e+00]])
        self.state= state_true
        # diff = state_true - ideal_state
        # print("Différence terme à terme (state_bis - self.state) :")
        # for i, (row_bis, row_state, row_diff) in enumerate(zip(state_true, self.state, diff)):
        #     print(f"Ligne {i}:")
        #     print(f"  state_bis : {row_bis}")
        #     print(f"  self.state: {row_state}")
        #     print(f"  Diff      : {row_diff}")
        self.previous_x = self.x
        
    def __call__(self,new_action=True):
        if new_action:
            
            self.state_list.append(self.state)
            print("State :",self.state)
            action = self.agent.select_action(self.state)
            self.past_action=action
            self.paraview_pos('paraview_export/') 
            
        else : 
            action = self.past_action
        action_global = coordinate_in_global_ref_3D(np.zeros(3), action,self.t, self.n, self.b)
        return action_global

    def history_state(self):
        if len(self.state_list)>0:
            path_save_fig = 'fig/'
            os.makedirs(path_save_fig, exist_ok=True)
            state_episode = np.array(self.state_list)
            n_steps = state_episode.shape[0]
            n_vec = state_episode.shape[1]
            states_reshaped = state_episode.reshape(n_steps, n_vec, 3)
            norms = np.linalg.norm(states_reshaped, axis=2)  # shape: (n_steps, n_vec)
            fig, axs = plt.subplots(2, 4, figsize=(15, 6))
            axs = axs.flatten()

            labels = [
                'Position local frame',
                'Velocity local frame',
                'Point 1 ahead',
                'Point 2 ahead',
                'Point 3 ahead',
                'Point 4 ahead',
                'Point 5 ahead',
                'Previous action'
            ]

            for i in range(n_vec):
                axs[i].hist(norms[:, i], bins=30, edgecolor='black')
                axs[i].set_title(labels[i])
                axs[i].set_xlabel("Magnitude")
                axs[i].set_ylabel("Frequency")
            plt.tight_layout()
            fig.suptitle(f"Distribution of State Vector Norms, distance : {self.d}", fontsize=16)
            plt.subplots_adjust(top=0.88)
            plt.savefig(os.path.join(path_save_fig, 'hist_norm.png'))
            plt.close(fig)
            
            ### Save state : 
            file_name = 'save_states/'
            os.makedirs(file_name, exist_ok=True)
            with open(os.path.join(file_name, "states.pkl"), "wb") as f:
                pickle.dump(self.state_list, f)
            

    def paraview_pos(self, output_save_path):
        os.makedirs(output_save_path, exist_ok=True)
        fig_path = 'fig/'
        os.makedirs(fig_path, exist_ok=True)

        # Save pos_wo_list as PolyData
        if self.pos_wo_list is not None and len(self.pos_wo_list) > 0:
            path_polydata = pv.PolyData(self.pos_wo_list)
            if len(self.pos_wo_list) > 1:
                lines = []
                for i in range(len(self.pos_wo_list) - 1):
                    lines.extend([2, i, i + 1])
                path_polydata.lines = np.array(lines)
                print(f"Nombre de segments créés: {len(self.pos_wo_list) - 1}")
            path_output = os.path.join(output_save_path, 'mean_positions_abf.vtp')
            path_polydata.save(path_output)
            print(f"Chemin sauvegardé: {path_output}")

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

        # Plot both histograms on the same figure (subplot)
        pos_list_arr = np.array(self.pos_wo_local_list) if self.pos_wo_local_list is not None and len(self.pos_wo_local_list) > 0 else None
        x_arr = np.array(self.x_list_local) if self.x_list_local is not None and len(self.x_list_local) > 0 else None
        pos_wo_arr = pos_list_arr
        if pos_wo_arr is not None or x_arr is not None:
            fig, axs = plt.subplots(1, 2, figsize=(12, 5))
            if pos_wo_arr is not None:
                norms_pos_wo = np.linalg.norm(pos_wo_arr, axis=1)
                axs[0].hist(norms_pos_wo, bins=30, edgecolor='black')
                axs[0].set_title('Normes des vecteurs pos_wo_list')
                axs[0].set_xlabel('Norme')
                axs[0].set_ylabel('Fréquence')
            else:
                axs[0].set_visible(False)
            if x_arr is not None:
                norms_x = np.linalg.norm(x_arr, axis=1)
                axs[1].hist(norms_x, bins=30, edgecolor='black')
                axs[1].set_title('Normes des vecteurs x_list')
                axs[1].set_xlabel('Norme')
                axs[1].set_ylabel('Fréquence')
            else:
                axs[1].set_visible(False)
            plt.tight_layout()
            plt.savefig(os.path.join(fig_path, 'hist_norm_poswo_xlist.png'))
            plt.close()

        # Plot trajectory in the yz-plane
        if x_arr is not None and x_arr.shape[1] >= 3:
            plt.figure(figsize=(8, 6))
            plt.plot(x_arr[:, 1], x_arr[:, 2], marker='o', linestyle='-', label='Trajectory (yz)')
            plt.plot(pos_wo_arr[:, 1], pos_wo_arr[:, 2], marker='o', linestyle='-', label='Trajectory (yz) filtered')
            plt.plot(x_arr[0,1],x_arr[0,2],marker='x',label ='Start')
            plt.plot(x_arr[-1,1],x_arr[-1,2],marker='x',label ='Last')
            plt.xlabel('y')
            plt.ylabel('z')
            plt.title('Trajectory in plane yz')
            plt.legend()
            plt.grid(True)
            plt.savefig(os.path.join(fig_path, 'trajectory_yz.png'))
            plt.close()

    def reward(self,x_target,dt):
        d = self.d
        rew_t = -self.C * dt
        rew_target = -np.linalg.norm(self.x - x_target) + np.linalg.norm(
            self.previous_x - x_target
        )
        rew_d = -self.beta * d
        rew = rew_t + rew_d + rew_target
        return rew_t, rew_d, rew_target, rew
    
def load_policy(policy_file: str, device_id: str, n_lookahead: int,previous_action:bool) -> callable:
    return PolicyEvaluator(policy_file, device_id, n_lookahead,previous_action)

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


