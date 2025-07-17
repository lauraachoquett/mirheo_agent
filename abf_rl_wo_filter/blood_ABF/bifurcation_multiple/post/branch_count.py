#!/usr/bin/env python

import numpy as np
import os
import pandas as pd
import re
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


def count_branching(basedir: str):
    params = load_parameters(os.path.join(basedir, 'parameters.pkl'))
    objstats_fname = os.path.join(basedir, 'obj_stats', 'ABF.csv')
    df = pd.read_csv(objstats_fname)

    ids = df['objId'].to_numpy()
    comx = df['comx'].to_numpy()
    comy = df['comy'].to_numpy()
    comz = df['comz'].to_numpy()
    time = df['time'].to_numpy()

    num_OBJs = len(np.unique(ids))
    num_steps = len(ids) // num_OBJs

    Lx, Ly, Lz = params.domain

    num_hi = 0
    num_lo = 0

    for i in range(num_OBJs):
        idx = np.argwhere(ids == i).flatten()
        x = comx[idx]
        y = comy[idx] - Ly/2
        z = comz[idx] - Lz/2
        t = time[idx]

        idx = np.argwhere(np.logical_and(x[:-1] < Lx/2, x[1:] > Lx/2)).flatten()

        for i in idx:
            if y[i] < 0:
                num_lo += 1
            else:
                num_hi += 1

    return num_lo, num_hi


def get_theta(basedir: str):
    rexf = '[-+]?\d*\.\d+|\d+'
    matches = re.findall(f"theta_({rexf})_", basedir)
    assert len(matches) == 1
    return float(matches[0])


def main(basedirs: list,
         out: str):

    thetas = []
    nums_lo = []
    nums_hi = []

    for basedir in basedirs:
        num_lo, num_hi = count_branching(basedir)
        theta = get_theta(basedir)

        nums_lo.append(num_lo)
        nums_hi.append(num_hi)
        thetas.append(theta)


    data = {'theta': thetas,
            'lower_branch': nums_lo,
            'upper_branch': nums_hi}

    if out is not None:
        df = pd.DataFrame(data)
        df.to_csv(out, index=False)
    else:
        print(data)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('basedirs', type=str, nargs='+', help="path of the simulation folder(s).")
    parser.add_argument('--out', type=str, default=None, help="output csv file.")
    args = parser.parse_args()

    main(args.basedirs, args.out)
