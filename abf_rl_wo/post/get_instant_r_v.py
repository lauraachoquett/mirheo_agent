#!/usr/bin/env python

import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import sys

here = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(here, ".."))

from cylinder.parameters import (CapilarryFlowParams,
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

def get_radial_traj(obj_stats_csv: str,
                    params_pkl: str):

    p = load_parameters(params_pkl)
    df = pd.read_csv(obj_stats_csv)

    obj_ids = df['objId'].to_numpy()
    nswimmers = len(np.unique(obj_ids))
    if nswimmers != 1:
        raise RuntimeError(f"got {nswimmers} swimmers in {obj_stats_csv}; expected 1.")

    t = df['time'].to_numpy()
    x = df['comx'].to_numpy()
    y = df['comy'].to_numpy()
    z = df['comz'].to_numpy()

    Lx = p.L
    Ly = 2 * (p.R + 2 * p.rc)
    Lz = 2 * (p.R + 2 * p.rc)

    x = unroll_periodic(x, Lx)

    y -= Ly/2
    z -= Lz/2
    r = np.sqrt(y**2 + z**2)

    v = np.diff(x) / np.diff(t)

    # non dimensionalize
    t *= p.omega
    r /= p.R
    v /= p.mean_vel

    v = np.append(v, v[-1]) # just to make it same length as others

    return t, r, v

def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('paths', type=str, nargs='+', help="path of the simulation folder(s).")
    parser.add_argument('--out', type=str, default=None, help="output csv file.")
    args = parser.parse_args(argv)

    all_t = []
    all_r = []
    all_v = []

    for basedir in args.paths:
        obj_stats_csv = os.path.join(basedir, 'obj_stats', 'ABF.csv')
        params_pkl = os.path.join(basedir, 'parameters.pkl')

        t, r, v = get_radial_traj(obj_stats_csv, params_pkl)
        all_t.append(t)
        all_r.append(r)
        all_v.append(v)

    niter = min([len(t) for t in all_t]) - 1

    data = {"tomega": all_t[0][:niter]}

    for i, r in enumerate(all_r):
        data[f"r{i}"] = r[:niter]

    for i, v in enumerate(all_v):
        data[f"v{i}"] = v[:niter]

    if args.out is not None:
        df = pd.DataFrame(data)
        df.to_csv(args.out, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
