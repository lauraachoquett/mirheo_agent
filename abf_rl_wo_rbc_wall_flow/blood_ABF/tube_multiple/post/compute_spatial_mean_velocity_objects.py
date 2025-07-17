#!/usr/bin/env python

import numpy as np
import os
import pandas as pd
import sys
from scipy.stats import binned_statistic

here = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(here, ".."))

from cylinder.parameters import (SimulationParams,
                        load_parameters)

def unroll_periodic(x: np.ndarray,
                    L: float):
    dx = np.diff(x)
    idx_pos = np.argwhere(dx > +L / 2).flatten()
    idx_neg = np.argwhere(dx < -L / 2).flatten()

    correction = np.zeros_like(x)
    correction[idx_neg] += L
    correction[idx_pos] -= L
    correction = np.cumsum(correction)

    x[1:] += correction[:-1]
    return x


def main(basedir: str,
         obj: str,
         out: str):
    params = load_parameters(os.path.join(basedir, 'parameters.pkl'))
    objstats_fname = os.path.join(basedir, 'stats', obj+'.csv')
    df = pd.read_csv(objstats_fname)

    ids = df['objId'].to_numpy()
    comx = df['comx'].to_numpy()
    comy = df['comy'].to_numpy()
    comz = df['comz'].to_numpy()

    num_OBJs = len(np.unique(ids))
    num_steps = len(ids) // num_OBJs

    volume_domain = params.L * np.pi * params.R**2

    if obj == 'ABF':
        volume_obj = params.mesh_ABF.volume
    else:
        volume_obj = params.mesh_ini.volume

    print(f"volume fraction: {num_OBJs * volume_obj / volume_domain} ({num_OBJs} {obj}s)")

    # use polar coordinates
    x = np.zeros((num_OBJs, num_steps))
    r = np.zeros((num_OBJs, num_steps))
    theta = np.zeros((num_OBJs, num_steps))

    cy = 2*params.rc + params.R
    cz = 2*params.rc + params.R

    for i in range(num_OBJs):
        idx = np.argwhere(ids == i).flatten()
        x[i,:] = unroll_periodic(comx[idx], params.L)

        y = comy[idx] - cy
        z = comz[idx] - cz
        r[i,:] = np.sqrt(y**2 + z**2)
        theta[i,:] = unroll_periodic(np.arctan2(y,z), 2 * np.pi)

    start_time = num_steps // 2
    x = x[:,start_time:]
    r = r[:,start_time:]
    theta = theta[:,start_time:]

    dx = np.diff(x, axis=1)
    dtheta = np.diff(theta, axis=1)
    r = (r[:,1:] + r[:,:-1]) / 2

    dt = 0.1 * 2 * np.pi / params.omega

    v_sim_to_mms = params.length_scale_.to('mm').magnitude / params.time_scale_.to('s').magnitude
    vx = dx/dt * v_sim_to_mms
    vt = -r*dtheta/dt * v_sim_to_mms
    r = r / params.R

    bins = np.linspace(0.0, 1.0, 10, endpoint=True)

    vx_median, re, _ = binned_statistic(r.flatten(), values=vx.flatten(), bins=bins, statistic='median')
    vt_median, re, _ = binned_statistic(r.flatten(), values=vt.flatten(), bins=bins, statistic='median')

    vx_mean, _, _ = binned_statistic(r.flatten(), values=vx.flatten(), bins=bins, statistic='mean')
    vt_mean, _, _ = binned_statistic(r.flatten(), values=vt.flatten(), bins=bins, statistic='mean')

    # 90% conf interval
    vx_95, _, _ = binned_statistic(r.flatten(), values=vx.flatten(), bins=bins, statistic=lambda x: np.quantile(x, q=0.95))
    vx_05, _, _ = binned_statistic(r.flatten(), values=vx.flatten(), bins=bins, statistic=lambda x: np.quantile(x, q=0.05))

    vt_95, _, _ = binned_statistic(r.flatten(), values=vt.flatten(), bins=bins, statistic=lambda x: np.quantile(x, q=0.95))
    vt_05, _, _ = binned_statistic(r.flatten(), values=vt.flatten(), bins=bins, statistic=lambda x: np.quantile(x, q=0.05))

    # 50% conf interval
    vx_75, _, _ = binned_statistic(r.flatten(), values=vx.flatten(), bins=bins, statistic=lambda x: np.quantile(x, q=0.75))
    vx_25, _, _ = binned_statistic(r.flatten(), values=vx.flatten(), bins=bins, statistic=lambda x: np.quantile(x, q=0.25))

    vt_75, _, _ = binned_statistic(r.flatten(), values=vt.flatten(), bins=bins, statistic=lambda x: np.quantile(x, q=0.75))
    vt_25, _, _ = binned_statistic(r.flatten(), values=vt.flatten(), bins=bins, statistic=lambda x: np.quantile(x, q=0.25))

    r = (re[1:] + re[:-1]) / 2

    df = pd.DataFrame({'r': r,
                       'vx_mean': vx_mean,
                       'vx_median': vx_median,
                       'vx_90_lo': vx_05,
                       'vx_90_hi': vx_95,
                       'vx_50_lo': vx_25,
                       'vx_50_hi': vx_75,
                       'vt_mean': vt_mean,
                       'vt_median': vt_median,
                       'vt_90_lo': vt_05,
                       'vt_90_hi': vt_95,
                       'vt_50_lo': vt_25,
                       'vt_50_hi': vt_75})
    if out is not None:
        df.to_csv(out, index=False)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('basedir', type=str, help="path of the simulation folder.")
    parser.add_argument('obj', type=str, choices=['ABF', 'RBC'], help="The kind of objects to measure.")
    parser.add_argument('--out', type=str, default=None, help="output csv file.")
    args = parser.parse_args()

    main(args.basedir, args.obj, args.out)
