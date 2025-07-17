s#!/usr/bin/env python

import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import sys

here = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(here, ".."))

from cylinder.parameters import (SimulationParams,
                        load_parameters)

def compute_mean_distance(obj_stats_csv: str,
                          params_pkl: str):

    p = load_parameters(params_pkl)
    df = pd.read_csv(obj_stats_csv)

    obj_ids = df['objId'].to_numpy()
    nswimmers = len(np.unique(obj_ids))
    if nswimmers != 1:
        raise RuntimeError(f"got {nswimmers} swimmers in {obj_stats_csv}; expected 1.")

    t = df['time'].to_numpy()
    y = df['comy'].to_numpy()
    z = df['comz'].to_numpy()

    Ly = 2 * (p.R + 2 * p.rc)
    Lz = 2 * (p.R + 2 * p.rc)

    y -= Ly/2
    z -= Lz/2
    r = np.sqrt(y**2 + z**2)
    rmean = np.mean(r)

    rdimless = rmean / p.R
    rdim = rmean * p.length_scale_
    return rdimless, rdim

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
    parser.add_argument('--against', type=str, default="a", help="Varying input parameter.")
    args = parser.parse_args(argv)

    param_name = args.against
    all_against_params = list()
    all_rdimless = list()
    all_rdim = list()

    for basedir in args.paths:
        obj_stats_csv = os.path.join(basedir, 'obj_stats', 'ABF.csv')
        params_pkl = os.path.join(basedir, 'parameters.pkl')

        rdimless, rdim = compute_mean_distance(obj_stats_csv, params_pkl)
        rdim = rdim.to('um').magnitude
        param_value = get_param(basedir, param_name)
        print(f"{param_name} = {param_value:.2e}, r = {rdimless:.4e} ({rdim} um)")
        all_rdimless.append(rdimless)
        all_rdim.append(rdim)
        all_against_params.append(param_value)

    if args.out is not None:
        df = pd.DataFrame({param_name : all_against_params,
                           "rdimless": all_rdimless,
                           "rdim": all_rdim})
        df.to_csv(args.out, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
