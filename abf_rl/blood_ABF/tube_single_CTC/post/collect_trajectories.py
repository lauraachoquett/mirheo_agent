#!/usr/bin/env python

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

def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('basedirs', type=str, nargs='+', help="Simulation directory(ies).")
    parser.add_argument('--out', type=str, required=True, help="The name of the output csv file.")
    args = parser.parse_args(argv)

    basedirs = args.basedirs

    all_times = []
    all_x_ABF = []
    all_x_CTC = []

    for i, basedir in enumerate(basedirs):
        ABF_csv = os.path.join(basedir, 'obj_stats', 'ABF.csv')
        CTC_csv = os.path.join(basedir, 'obj_stats', 'CTC.csv')
        params = os.path.join(basedir, 'parameters.pkl')

        p = load_parameters(params)
        df_ABF = pd.read_csv(ABF_csv)
        df_CTC = pd.read_csv(CTC_csv)

        time = df_ABF['time'].to_numpy()
        x_ABF = df_ABF['comx'].to_numpy()
        x_CTC = df_CTC['comx'].to_numpy()

        x_ABF = unroll_periodic(x_ABF, p.L) / p.ABF_length
        x_CTC = unroll_periodic(x_CTC, p.L) / p.ABF_length
        time *= p.omega

        all_times.append(time)
        all_x_ABF.append(x_ABF)
        all_x_CTC.append(x_CTC)

    end = min([len(time) for time in all_times])

    data = dict()
    data['time'] = all_times[0][:end]

    for i in range(len(basedirs)):
        data[f"x_ABF_{i}"] = all_x_ABF[i][:end]
        data[f"x_CTC_{i}"] = all_x_CTC[i][:end]

    out_csv = args.out

    if out_csv is not None:
        df = pd.DataFrame(data)
        df.to_csv(out_csv, index=False)

if __name__ == '__main__':
    main(sys.argv[1:])
