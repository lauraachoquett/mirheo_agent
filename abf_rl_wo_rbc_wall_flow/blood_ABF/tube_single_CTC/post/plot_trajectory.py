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

def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('basedirs', type=str, nargs='+', help="Simulation directory(ies).")
    parser.add_argument('--out', type=str, default='GUI', help="The name of the output file.")
    parser.add_argument('--show-diff', action='store_true', default=False, help="Show difference in x positions of CTC and ABF.")
    args = parser.parse_args(argv)

    basedirs = args.basedirs

    fig, ax = plt.subplots(figsize=[6.8, 6.8])

    for basedir in basedirs:
        ABF_csv = os.path.join(basedir, 'obj_stats', 'ABF.csv')
        CTC_csv = os.path.join(basedir, 'obj_stats', 'CTC.csv')
        params = os.path.join(basedir, 'parameters.pkl')

        p = load_parameters(params)
        df_ABF = pd.read_csv(ABF_csv)
        df_CTC = pd.read_csv(CTC_csv)

        time = df_ABF['time'].to_numpy()
        x_ABF = df_ABF['comx'].to_numpy()
        x_CTC = df_CTC['comx'].to_numpy()

        x_ABF = unroll_periodic(x_ABF, p.L) / p.L
        x_CTC = unroll_periodic(x_CTC, p.L) / p.L
        time *= p.omega

        if args.show_diff:
            ax.plot(time, x_CTC - x_ABF)
        else:
            xlo = x_ABF - p.ABF_length/2 / p.L
            xhi = x_ABF + p.ABF_length/2 / p.L
            ax.fill_between(time, xlo, xhi, alpha=0.5, label="ABF")
            #ax.plot(time, x_ABF, label="ABF")

            xlo = x_CTC - p.CTC_diameter/2 / p.L
            xhi = x_CTC + p.CTC_diameter/2 / p.L
            ax.fill_between(time, xlo, xhi, alpha=0.5, label="CTC")
            #ax.plot(time, x_CTC, label="CTC")
            ax.legend(loc='lower right')

    ax.set_xlim(0,np.max(time))
    ax.set_ylim(0,)

    ax.set_xlabel(f'$\omega t$')
    ax.set_ylabel(f'$x/L$')

    if args.out == 'GUI':
        plt.show()
    else:
        plt.savefig(args.out, transparent=True)

if __name__ == '__main__':
    main(sys.argv[1:])
