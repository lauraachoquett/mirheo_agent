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

    x = np.zeros((num_OBJs, num_steps))
    r = np.zeros((num_OBJs, num_steps))

    cy = 2*params.rc + params.R
    cz = 2*params.rc + params.R

    for i in range(num_OBJs):
        idx = np.argwhere(ids == i).flatten()
        x[i,:] = unroll_periodic(comx[idx], params.L)

        y = comy[idx] - cy
        z = comz[idx] - cz
        r[i,:] = np.sqrt(y**2 + z**2)

    r = r / params.R
    bins = np.linspace(0.0, 1.0, 10, endpoint=True)
    half_time = num_steps // 2

    # first half of simulation time
    r_ = r[:,:half_time]
    rho1, re = np.histogram(r_.flatten(), bins=bins)
    volumes = np.pi * params.R**2 * (re[1:]**2 - re[:-1]**2) * params.L * params.length_scale_.to('um')**3
    rho1 = rho1 / (volumes * r_.shape[1])


    # second half of simulation time
    r_ = r[:,half_time:]
    rho2, re = np.histogram(r_.flatten(), bins=bins)
    volumes = np.pi * params.R**2 * (re[1:]**2 - re[:-1]**2) * params.L * params.length_scale_.to('um')**3
    rho2 = rho2 / (volumes * r_.shape[1])

    r = (re[1:] + re[:-1]) / 2

    df = pd.DataFrame({'r': r,
                       'rho1': rho1,
                       'rho2': rho2})

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
