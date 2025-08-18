#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

def plot(fname: str):
    df = pd.read_csv(fname)

    fig, axes = plt.subplots(ncols=2, nrows=1,
                             figsize=(6.4*2, 4.8))

    ax = axes[0]
    ax.plot(df['time'], df['fx'], '-o')
    ax.plot(df['time'], df['fy'], '-o')
    ax.plot(df['time'], df['fz'], '-o')

    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$f$')
    ax.legend()

    ax = axes[1]
    ax.plot(df['time'], df['Tx'], '-o')
    ax.plot(df['time'], df['Ty'], '-o')
    ax.plot(df['time'], df['Tz'], '-o')

    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$f$')
    ax.legend()

    plt.show()

def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('stats', type=str, help="stats file.")
    args = parser.parse_args(argv)

    plot(args.stats)

if __name__ == '__main__':
    main(sys.argv[1:])
