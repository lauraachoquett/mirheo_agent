#!/usr/bin/env python

import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import sys


def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('csv', type=str, help="Csv file that contains the wobbling angle and Hematocrit for several simulations.")
    parser.add_argument('--out', type=str, default='GUI', help="The name of the output file.")
    args = parser.parse_args(argv)

    df = pd.read_csv(args.csv)

    all_Ht = df['Ht'].to_numpy()
    all_thetas = df['theta'].to_numpy()

    Ht = sorted(np.unique(all_Ht))
    theta = list()
    dtheta = list()

    for h in Ht:
        idx = np.argwhere(all_Ht == h).flatten()
        thetas = all_thetas[idx]
        theta.append(np.mean(thetas))
        dtheta.append(np.std(thetas) / np.sqrt(len(idx)))

    Ht = np.array(Ht)
    theta = np.array(theta)
    dtheta = np.array(dtheta)

    fig, ax = plt.subplots()
    ax.plot(all_Ht, all_thetas, '.')
    ax.errorbar(Ht, theta, yerr=dtheta, fmt='ko', capsize=3)
    ax.set_xlabel(r'$Ht$')
    ax.set_ylabel(r'$\theta$')

    ax.set_ylim(0,)

    if args.out == 'GUI':
        plt.show()
    else:
        plt.savefig(args.out, transparent=True)

if __name__ == '__main__':
    main(sys.argv[1:])
