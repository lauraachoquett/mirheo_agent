#!/usr/bin/env python3

import json
import numpy as np
from TD3 import TD3
from utils import coordinate_in_path_ref_3D,coordinate_in_global_ref_3D
import torch
import matplotlib.pyplot as plt
import os 
import pyvista as pv


from scipy.signal import butter, filtfilt

class RobustPositionFilter:
    def __init__(self, omega_precession, dt, order=2):
        f_precession = omega_precession / (2 * np.pi)
        self.fc = f_precession / 5
        self.dt = dt
        self.order = order
        self.omega_precession = omega_precession
        cycles_needed = 6
        self.buffer_size_freq = int(cycles_needed * 2*np.pi / (self.omega_precession * self.dt))

        nyquist = 0.5 / dt
        normalized_fc = self.fc / nyquist
        self.b, self.a = butter(order, normalized_fc, btype='low')
        
        self.x_buffer = []
        self.y_buffer = []
        self.z_buffer = []
        
    def filter_position(self, new_position):
        self.x_buffer.append(new_position[0])
        self.y_buffer.append(new_position[1])
        self.z_buffer.append(new_position[2])
        
        buffer_size = max(100,self.buffer_size_freq)

        if len(self.x_buffer) > buffer_size:
            self.x_buffer.pop(0)
            self.y_buffer.pop(0)
            self.z_buffer.pop(0)
        
        if len(self.x_buffer) >= 2*self.order + 6:
            x_filt = filtfilt(self.b, self.a, self.x_buffer)[-1]
            y_filt = filtfilt(self.b, self.a, self.y_buffer)[-1]
            z_filt = filtfilt(self.b, self.a, self.z_buffer)[-1]
            return np.array([x_filt, y_filt, z_filt])
        else:
            return new_position
        
class PolicyEvaluator:
    def __init__(self, policy_file: str, device_id: str, n_lookahead: int,omega:float,dt:float):
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
        self.typical_length_rom = 2
        self.length_scale = self.length_cylinder / self.typical_length_rom
        
        self.time_rom = 3.333e-02 / 5
        self.time_sim = 0.00163
        self.time_scale = self.time_sim / self.time_rom
        self.velocity_scale = 1.0 / 5e-3 
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
        self.x_list_local=[]

        self.x_list_filter =[]
        self.x_list_local_filter = []
        self.position_filter = RobustPositionFilter(omega_precession=omega, dt=dt)
    
    
    def init_state(self, pos, path, T, N, B):
        self.x_list.append(pos)
        self.x_list_filter.append(pos)
        self.x = pos
        self.previous_x = self.x
        self.past_action = np.array(T[0])
        path = path
        d, id_cp = min_dist_closest_point(self.x, path)
        self.d=d
        self.p_cp = path[id_cp]
        self.t = T[id_cp]
        self.n = N[id_cp]
        self.b = B[id_cp]
        
    def compute_state(self, pos, path, T, N, B, dt):
        filtered_position = self.position_filter.filter_position(pos)
        
        self.x_list.append(pos)
        self.x_list_filter.append(filtered_position)
        
        self.x = filtered_position
        d, id_cp = min_dist_closest_point(self.x, path)
        self.d=d
        path_len = len(path)
        self.p_cp = path[id_cp]
        self.t = T[id_cp]
        self.n = N[id_cp]
        self.b = B[id_cp]

        s = coordinate_in_path_ref_3D(self.p_cp, self.x, self.t, self.n, self.b)
        self.x_list_local_filter.append(s)
        self.x_list_local.append(coordinate_in_path_ref_3D(self.p_cp, pos, self.t, self.n, self.b))
        
        s_scaled = s/self.length_scale
        result = [s_scaled.reshape(1, 3)]
        
        s_previous = coordinate_in_path_ref_3D(self.p_cp, self.previous_x, self.t, self.n, self.b)
        v_local_path = (s - s_previous) / dt
        v_norm = np.linalg.norm(v_local_path) * self.velocity_scale
        print("V_norm : ",v_norm)
        v_local_path = v_local_path *  self.velocity_scale
        if v_norm > 5 and v_norm != 0:
            v_local_path = v_local_path / v_norm * 5
        result.append(v_local_path .reshape(1, 3))

        if self.n_lookahead > 0:
            lookahead = []
            for i in range(1, self.n_lookahead + 1):
                idx = min(id_cp + i, path_len - 1)
                next_p = path[idx]
                next_p_local = coordinate_in_path_ref_3D(self.p_cp, next_p, self.t, self.n, self.b).reshape(1, 3)
                lookahead.append(next_p_local/self.length_scale)
            result.append(np.concatenate(lookahead, axis=0))
        result.append(self.past_action.reshape(1, 3))
        self.state = np.concatenate(result, axis=0)

        self.previous_x=self.x
        
    def __call__(self,new_action=True):
        if new_action:
            self.state_list.append(self.state)
            action = self.agent.select_action(self.state)
            self.past_action=action
            print("Past action : ",self.past_action)
            print("Norm past action :", np.linalg.norm(self.past_action))
            self.paraview_pos('paraview_export/') 
            
        else : 
            action = self.past_action
        action_global = coordinate_in_global_ref_3D(np.zeros(3), action,self.t, self.n, self.b)
        return action_global

    def history_state(self):
        path_save_fig = 'fig/'
        os.makedirs(path_save_fig, exist_ok=True)
        state_episode = np.array(self.state_list)
        n_steps = state_episode.shape[0]
        n_vec = 8
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

    def paraview_pos(self, output_save_path):
        os.makedirs(output_save_path, exist_ok=True)
        fig_path = 'fig/'
        os.makedirs(fig_path, exist_ok=True)

        # Save pos_wo_list as PolyData
        if self.x_list_filter is not None and len(self.x_list_filter) > 0:
            path_polydata = pv.PolyData(self.x_list_filter)
            if len(self.x_list_filter) > 1:
                lines = []
                for i in range(len(self.x_list_filter) - 1):
                    lines.extend([2, i, i + 1])
                path_polydata.lines = np.array(lines)
                print(f"Nombre de segments créés: {len(self.x_list_filter) - 1}")
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
        pos_list_arr = np.array(self.x_list_local_filter) if self.x_list_local_filter is not None and len(self.x_list_local_filter) > 0 else None
        x_arr = np.array(self.x_list_local) if self.x_list_local is not None and len(self.x_list_local) > 0 else None
        pos_wo_arr = pos_list_arr
        if pos_wo_arr is not None or x_arr is not None:
            fig, axs = plt.subplots(1, 2, figsize=(12, 5))
            if pos_wo_arr is not None:
                norms_pos_wo = np.linalg.norm(pos_wo_arr, axis=1)
                axs[0].hist(norms_pos_wo, bins=30, edgecolor='black')
                axs[0].set_title('Magnitude local filetered position ')
                axs[0].set_xlabel('Magnitude')
                axs[0].set_ylabel('Frequency')
            else:
                axs[0].set_visible(False)
            if x_arr is not None:
                norms_x = np.linalg.norm(x_arr, axis=1)
                axs[1].hist(norms_x, bins=30, edgecolor='black')
                axs[1].set_title('Magnitude local true position')
                axs[1].set_xlabel('Magnitude')
                axs[1].set_ylabel('Frequency')
            else:
                axs[1].set_visible(False)
            plt.tight_layout()
            plt.savefig(os.path.join(fig_path, 'hist_norm_poswo_xlist.png'))
            plt.close()

        # Plot trajectory in the yz-plane
        if x_arr is not None and x_arr.shape[1] >= 3:
            plt.figure(figsize=(10, 8))
            # Couleurs cohérentes pour chaque trajectoire
            color_brute = '#1f77b4'
            color_filtre = '#ff7f0e'

            # Trajectoire brute
            plt.plot(
            x_arr[:, 1], x_arr[:, 2],
            marker='o', linestyle='-', color=color_brute, alpha=0.7, linewidth=2, label='Trajectoire brute (yz)'
            )
            # Début et fin brute (même couleur que la courbe)
            plt.scatter(x_arr[0, 1], x_arr[0, 2], marker='*', color=color_brute, s=180, edgecolor='black', label='Départ brute')
            plt.scatter(x_arr[-1, 1], x_arr[-1, 2], marker='X', color=color_brute, s=180, edgecolor='black', label='Arrivée brute')

            # Trajectoire filtrée
            if pos_wo_arr is not None and pos_wo_arr.shape[1] >= 3:
                plt.plot(
                    pos_wo_arr[:, 1], pos_wo_arr[:, 2],
                    marker='s', linestyle='-', color=color_filtre, alpha=0.8, linewidth=2, label='Trajectoire filtrée (yz)'
                )
                # Début et fin filtrée (même couleur que la courbe)
                plt.scatter(pos_wo_arr[0, 1], pos_wo_arr[0, 2], marker='*', color=color_filtre, s=160, edgecolor='black', label='Départ filtrée')
                plt.scatter(pos_wo_arr[-1, 1], pos_wo_arr[-1, 2], marker='X', color=color_filtre, s=160, edgecolor='black', label='Arrivée filtrée')

            # Axes et titre
            plt.xlabel('y', fontsize=14)
            plt.ylabel('z', fontsize=14)
            plt.title('Trajectoire dans le plan (y, z)', fontsize=16)
            plt.legend(fontsize=12, loc='best', frameon=True)
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.tight_layout()
            plt.savefig(os.path.join(fig_path, 'trajectory_yz.png'))
            plt.close()

        # Nouvelle figure : évolution de x en fonction du temps
        if x_arr is not None and x_arr.shape[1] >= 1:
            plt.figure(figsize=(10, 6))
            t = np.arange(x_arr.shape[0])
            plt.plot(t, x_arr[:, 0], label='x brute', color=color_brute, marker='o', alpha=0.7)
            if pos_wo_arr is not None and pos_wo_arr.shape[1] >= 1:
                plt.plot(t, pos_wo_arr[:, 0], label='x filtrée', color=color_filtre, marker='s', alpha=0.8)
            plt.xlabel('Temps (itération)', fontsize=14)
            plt.ylabel('x', fontsize=14)
            plt.title('Évolution de x en fonction du temps', fontsize=16)
            plt.legend(fontsize=12, loc='best', frameon=True)
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.tight_layout()
            plt.savefig(os.path.join(fig_path, 'x_vs_time.png'))
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
    
def load_policy(policy_file: str, device_id: str, n_lookahead: int,omega:float,dt:float) -> callable:
    return PolicyEvaluator(policy_file, device_id, n_lookahead,omega,dt)


def min_dist_closest_point(x, path):
    path = np.asarray(path)
    x = np.asarray(x)
    dists = np.linalg.norm(path - x, axis=-1)
    idx = np.argmin(dists)
    return dists[idx], idx


