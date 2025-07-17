#!/usr/bin/env python

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

def read_x(obj_stats_csv: str,
                  params_pkl: str):

    p = load_parameters(params_pkl)
    df = pd.read_csv(obj_stats_csv)

    obj_ids = df['objId'].to_numpy()
    nswimmers = len(np.unique(obj_ids))
    assert nswimmers == 1

    t = df['time'].to_numpy()
    x = df['comx'].to_numpy()

    x = unroll_periodic(x, p.L)

    x /= p.abf_length
    t *= p.omega

    return t, x

def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('paths', type=str, nargs='+', help="path of the simulation folder(s).")
    parser.add_argument('--out', type=str, default=None, help="if set, name of the csv file containing x/l vs t*omega")
    args = parser.parse_args(argv)

    all_t = list()
    all_x = list()

    for basedir in args.paths:
        obj_stats_csv = os.path.join(basedir, 'obj_stats', 'abf.csv')
        params_pkl = os.path.join(basedir, 'parameters.pkl')

        t, x = read_x(obj_stats_csv, params_pkl)

        all_t.append(t)
        all_x.append(x)

    n = min([len(t) for t in all_t])

    t = all_t[0][:n]
    all_x = [x[:n] for x in all_x]
    all_x = np.array(all_x)

    xmean = np.mean(all_x, axis=0)
    p = np.polyfit(t, xmean, deg=1)
    vx = p[0]
    print(f"Swimming speed: {vx} [body_length omega]")

    msd = np.mean((all_x - xmean[np.newaxis,:])**2, axis=0)
    p = np.polyfit(t, msd, deg=1)
    D = p[0] / 2
    print(f"Diffusion coeff: {D} [body_length**2 omega]")

    if args.out is not None:
        import pandas as pd
        data = {"t": t}
        for i, x in enumerate(all_x):
            data[f"x{i}"] = x-x[0]
        df = pd.DataFrame(data)
        df.to_csv(args.out, index=False)

    if False:
        msd_ = np.polyval(p, t)
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(t, msd, '-')
        ax.plot(t, msd_, '--k')
        plt.show()

if __name__ == '__main__':
    main(sys.argv[1:])
