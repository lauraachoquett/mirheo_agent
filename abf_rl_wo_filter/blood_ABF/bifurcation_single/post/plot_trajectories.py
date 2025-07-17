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
    parser.add_argument('paths', type=str, nargs='+', help="path of the simulation folder(s).")
    args = parser.parse_args(argv)

    fig, ax = plt.subplots()

    for basedir in args.paths:
        obj_stats_csv = os.path.join(basedir, 'obj_stats', 'ABF.csv')
        params_pkl = os.path.join(basedir, 'parameters.pkl')
        p = load_parameters(params_pkl)

        df = pd.read_csv(obj_stats_csv)
        x = df['comx'].to_numpy()
        y = df['comy'].to_numpy()

        x = unroll_periodic(x, L=p.domain[0])

        ax.plot(x, y, '-')

    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    plt.show()


if __name__ == '__main__':
    main(sys.argv[1:])
