#!/usr/bin/env python

import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import sys


def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('csv', type=str, help="Csv file that contains the mean velocity and Hematocrit for several simulations.")
    parser.add_argument('--out', type=str, default='GUI', help="The name of the output file.")
    args = parser.parse_args(argv)

    df = pd.read_csv(args.csv)

    all_a = df['a'].to_numpy()
    all_v = df['vdimless'].to_numpy()

    a = sorted(np.unique(all_a))
    v = list()
    dv = list()

    for alpha in a:
        idx = np.argwhere(all_a == alpha).flatten()
        vs = all_v[idx]
        v.append(np.mean(vs))
        dv.append(np.std(vs) / np.sqrt(len(idx)))

    a = np.array(a)
    v = np.array(v)
    dv = np.array(dv)

    fig, ax = plt.subplots()
    ax.plot(all_a, all_v, '.')
    ax.errorbar(a, v, yerr=dv, fmt='ko', capsize=3)
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$V / l \Omega$')

    ax.set_xscale('log')
    #ax.set_ylim(0,)

    if args.out == 'GUI':
        plt.show()
    else:
        plt.savefig(args.out, transparent=True)

if __name__ == '__main__':
    main(sys.argv[1:])
