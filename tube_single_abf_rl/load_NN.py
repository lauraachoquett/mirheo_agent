#!/usr/bin/env python3

import json
import numpy as np
from TD3 import TD3
from utils import coordinate_in_path_ref_3D,coordinate_in_global_ref_3D
import torch
import matplotlib.pyplot as plt
import os 
import pyvista as pv

class PolicyEvaluator:
    def __init__(self, policy_file: str, device_id: str, n_lookahead: int):
        self.device = torch.device(f"cuda:{device_id}")
        self.n_lookahead = n_lookahead

        state_dim = 3 * (1 + 1 + 1 + n_lookahead)
        action_dim = 3
        max_action = np.inf
        self.agent = TD3(state_dim, action_dim, max_action)
        self.agent.load(policy_file, device=self.device)

        self.radius_biggest_ca_p = 0.269
        self.radius_biggest_ca_grid = 6.43
        self.radius_biggest_phys = 20  # micrometer

        self.length_cylinder = 150
        self.typical_length_rom = 5
        self.length_scale = self.length_cylinder / self.typical_length_rom
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
        self.x_list=[]
        self.time_rom = 3.333e-02 / 5
        self.time_sim = 0.00163
        self.time_scale = self.time_sim / self.time_rom
        
    def compute_state(self, x, path, T, N, B, dt):
        self.x_list.append(x)
        self.x = x/self.length_scale 
        
        if self.previous_x is None:
            self.previous_x = self.x
        if self.past_action is None : 
            self.past_action = np.array(T[0])

        path = path/self.length_scale
        d, id_cp = min_dist_closest_point(self.x, path)
        self.d=d*self.length_scale 
        path_len = len(path)
        p_cp = path[id_cp]
        self.t = T[id_cp]
        self.n = N[id_cp]
        self.b = B[id_cp]

        s = coordinate_in_path_ref_3D(p_cp, self.x, self.t, self.n, self.b)
        result = [s.reshape(1, 3)]
        
        s_previous = coordinate_in_path_ref_3D(p_cp, self.previous_x, self.t, self.n, self.b)
        v_local_path = (s - s_previous) / dt
        v_norm = np.linalg.norm(v_local_path)
        if v_norm > 5 and v_norm != 0:
            v_local_path = v_local_path / v_norm * 5
        result.append(v_local_path.reshape(1, 3))

        if self.n_lookahead > 0:
            lookahead = []
            for i in range(1, self.n_lookahead + 1):
                idx = min(id_cp + i, path_len - 1)
                next_p = path[idx]
                next_p_local = coordinate_in_path_ref_3D(p_cp, next_p, self.t, self.n, self.b).reshape(1, 3)
                lookahead.append(next_p_local)
            result.append(np.concatenate(lookahead, axis=0))
        result.append(self.past_action.reshape(1, 3))
        self.state = np.concatenate(result, axis=0)

        self.previous_x=self.x
        
    def __call__(self,new_action=True):
        if new_action:
            self.state_list.append(self.state)
            action = self.agent.select_action(self.state)
            self.past_action=action
        else : 
            action = self.past_action
        # print("Action :",action)
        action_global = coordinate_in_global_ref_3D(np.zeros(3), action,self.t, self.n, self.b)
        # print("Action global ref :",action_global)
        return action_global

    def history_state(self):
        path_save_fig = 'fig/'
        os.makedirs(path_save_fig, exist_ok=True)
        state_episode = np.array(self.state_list)
        # Check the shape to avoid reshape errors
        n_steps = state_episode.shape[0]
        n_vec = 8
        states_reshaped = state_episode.reshape(n_steps, n_vec, 3)
        norms = np.linalg.norm(states_reshaped, axis=2)  # shape: (n_steps, n_vec)
        fig, axs = plt.subplots(2, 4, figsize=(15, 6))
        axs = axs.flatten()

        for i in range(n_vec):
            axs[i].hist(norms[:, i], bins=30, edgecolor='black')
            axs[i].set_title(f"Vector {i}")
            axs[i].set_xlabel("Magnitude")
            axs[i].set_ylabel("Frequency")
        plt.tight_layout()
        fig.suptitle(f"Distribution of State Vector Norms, distance : {self.d}", fontsize=16)  # Titre global
        plt.subplots_adjust(top=0.88)  # Ajuste l'espace pour le titre global
        plt.savefig(os.path.join(path_save_fig, 'hist_norm.png'))
        plt.close(fig)

    def paraview_pos(self,output_save_path):
        os.makedirs(output_save_path,exist_ok=True)
        if self.x_list is not None and len(self.x_list) > 0:
            path_polydata = pv.PolyData(self.x_list)
            
            if len(self.x_list) > 1:
                lines = []
                for i in range(len(self.x_list) - 1):
                    lines.extend([2, i, i + 1])  
                path_polydata.lines = np.array(lines)
                print(f"Nombre de segments créés: {len(self.x_list) - 1}")
            
            path_output = os.path.join(output_save_path,'positions_abf.vtp')
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
    
def load_policy(policy_file: str, device_id: str, n_lookahead: int) -> callable:
    return PolicyEvaluator(policy_file, device_id, n_lookahead)


def min_dist_closest_point(x, path):
    path = np.asarray(path)
    x = np.asarray(x)
    dists = np.linalg.norm(path - x, axis=-1)
    idx = np.argmin(dists)
    return dists[idx], idx


