#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot(filenames: list):
    fig, ax = plt.subplots()

    for filename in filenames:
        df = pd.read_csv(filename)
        r = df['r'].to_numpy()
        vx = df['vx_mean'].to_numpy()
        ax.plot(r, vx)

    ax.set_xlabel(r"$r/R$")
    ax.set_ylabel(r"$v_x$ [$\mu$m/s]")
    plt.show()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('csv', type=str, nargs='+', help="The data to plot.")
    args = parser.parse_args()

    plot(filenames=args.csv)
