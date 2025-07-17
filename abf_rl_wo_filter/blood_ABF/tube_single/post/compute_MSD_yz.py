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

def compute_MSD_yz(obj_stats_csv: str,
                   params_pkl: str):

    p = load_parameters(params_pkl)
    df = pd.read_csv(obj_stats_csv)

    obj_ids = df['objId'].to_numpy()
    nswimmers = len(np.unique(obj_ids))
    assert nswimmers == 1

    t = df['time'].to_numpy()
    y = df['comy'].to_numpy()
    z = df['comz'].to_numpy()


    yc = p.R + 2 * p.rc
    zc = p.R + 2 * p.rc
    dist_sq = (y-yc)**2 + (z-zc)**2

    return t, dist_sq


def get_Ht(basedir: str):
    import re
    rexf = '[-+]?\d*\.\d+|\d+'
    matches = re.findall(f"Ht_({rexf})", basedir)
    assert len(matches) == 1
    return float(matches[0])

def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('paths', type=str, nargs='+', help="path of the simulation folder(s).")
    parser.add_argument('--out', type=str, default="results/Ht_v.csv", help="output csv file.")
    args = parser.parse_args(argv)

    all_t = list()
    all_d2 = list()

    for basedir in args.paths:
        obj_stats_csv = os.path.join(basedir, 'obj_stats', 'ABF.csv')
        params_pkl = os.path.join(basedir, 'parameters.pkl')

        t, d2 = compute_MSD_yz(obj_stats_csv, params_pkl)

        all_t.append(t)
        all_d2.append(d2)

    n = min([len(t) for t in all_t])

    t = all_t[0][:n]
    msd = np.zeros_like(t)

    for d2 in all_d2:
        msd += d2[:n]

    msd /= len(all_d2)

    fig, ax = plt.subplots()
    ax.plot(t, np.sqrt(msd))
    plt.show()

if __name__ == '__main__':
    main(sys.argv[1:])
