#!/usr/bin/env python

import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import sys

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

def compute_velocity(obj_stats_csv: str,
                     params_pkl: str):

    p = load_parameters(params_pkl)
    df = pd.read_csv(obj_stats_csv)

    obj_ids = df['objId'].to_numpy()
    nswimmers = len(np.unique(obj_ids))
    if nswimmers != 1:
        raise RuntimeError(f"got {nswimmers} swimmers in {obj_stats_csv}; expected 1.")

    t = df['time'].to_numpy()
    x = df['comx'].to_numpy()
    x = unroll_periodic(x, p.L)

    #vmean = np.polyfit(t, x, deg=1)[0]
    vmean = (x[-1] - x[0]) / (t[-1] - t[0])

    vdimless = vmean / (p.abf_length * p.omega)
    vdim = vmean * p.length_scale_ / p.time_scale_
    return vdimless, vdim

def get_param(basedir: str, param_name: str):
    import re
    rexf = '[-+]?\d*\.\d+|\d+'
    matches = re.findall(f"{param_name}_({rexf})", basedir)
    assert len(matches) == 1
    return float(matches[0])

def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('paths', type=str, nargs='+', help="path of the simulation folder(s).")
    parser.add_argument('--out', type=str, default=None, help="output csv file.")
    parser.add_argument('--against', type=str, default="Ht", help="Varying input parameter.")
    args = parser.parse_args(argv)

    param_name = args.against
    all_against_params = list()
    all_vdimless = list()
    all_vdim = list()

    for basedir in args.paths:
        obj_stats_csv = os.path.join(basedir, 'obj_stats', 'abf.csv')
        params_pkl = os.path.join(basedir, 'parameters.pkl')

        vdimless, vdim = compute_velocity(obj_stats_csv, params_pkl)
        vdim = vdim.to('um/s').magnitude
        param_value = get_param(basedir, param_name)
        print(f"{param_name} = {param_value:.2e}, v = {vdimless:.4e} ({vdim} um/s)")
        all_vdimless.append(vdimless)
        all_vdim.append(vdim)
        all_against_params.append(param_value)

    if args.out is not None:
        df = pd.DataFrame({param_name : all_against_params,
                           "vdimless": all_vdimless,
                           "vdim": all_vdim})
        df.to_csv(args.out, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
