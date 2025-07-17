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

    all_Ht = df['Ht'].to_numpy()
    all_v = df['vdimless'].to_numpy()

    Ht = sorted(np.unique(all_Ht))
    v = list()
    dv = list()

    for h in Ht:
        idx = np.argwhere(all_Ht == h).flatten()
        vs = all_v[idx]
        v.append(np.mean(vs))
        dv.append(np.std(vs) / np.sqrt(len(idx)))

    Ht = np.array(Ht)
    v = np.array(v)
    dv = np.array(dv)

    fig, ax = plt.subplots()
    ax.plot(all_Ht, all_v, '.')
    ax.errorbar(Ht, v, yerr=dv, fmt='ko', capsize=3)
    ax.set_xlabel(r'$Ht$')
    ax.set_ylabel(r'$V / l \Omega$')

    ax.set_ylim(0,)

    if args.out == 'GUI':
        plt.show()
    else:
        plt.savefig(args.out, transparent=True)

if __name__ == '__main__':
    main(sys.argv[1:])
