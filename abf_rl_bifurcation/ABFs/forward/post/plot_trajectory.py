#!/usr/bin/env python

import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import sys

here = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(here, ".."))

from cylinder.parameters import (Parameters,
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
    parser.add_argument('basedir', type=str, help="Simulation directory.")
    parser.add_argument('--out', type=str, default='GUI', help="The name of the output file.")
    args = parser.parse_args(argv)

    basedir = args.basedir
    csv    = os.path.join(basedir, 'obj_stats', 'ABF.csv')
    params = os.path.join(basedir, 'parameters.pkl')

    p  = load_parameters(params)
    df = pd.read_csv(csv)

    t = df['time'].to_numpy()
    x = df['comx'].to_numpy()

    t *= p.omega
    x = unroll_periodic(x, p.domain[0]) / p.ABF_L

    fig, ax = plt.subplots()
    ax.plot(t, x)
    ax.set_xlabel(f'$\omega t$')
    ax.set_ylabel(f'$x/l$')

    if args.out == 'GUI':
        plt.show()
    else:
        plt.savefig(args.out, transparent=True)

if __name__ == '__main__':
    main(sys.argv[1:])
