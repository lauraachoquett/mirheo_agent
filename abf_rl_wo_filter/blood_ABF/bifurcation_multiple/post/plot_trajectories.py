#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
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




def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=str, help="path of the simulation folder.")
    args = parser.parse_args(argv)

    fig, ax = plt.subplots()

    basedir = args.path
    obj_stats_csv = os.path.join(basedir, 'obj_stats', 'ABF.csv')
    params_pkl = os.path.join(basedir, 'parameters.pkl')
    p = load_parameters(params_pkl)

    df = pd.read_csv(obj_stats_csv)
    ids = df['objId'].to_numpy()
    comx = df['comx'].to_numpy()
    comy = df['comy'].to_numpy()

    num_OBJs = len(np.unique(ids))
    num_steps = len(ids) // num_OBJs

    for i in range(num_OBJs):
        idx = np.argwhere(ids == i).flatten()
        x = unroll_periodic(comx[idx], L=p.domain[0])
        y = comy[idx]
        ax.plot(x, y, '-')

    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    plt.show()


if __name__ == '__main__':
    main(sys.argv[1:])
