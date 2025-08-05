#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys


def main(argv: list):
    parser = argparse.ArgumentParser()
    parser.add_argument('sim', type=str, help="Simulation results.")
    args = parser.parse_args(argv)

    fig, ax = plt.subplots()

    df = pd.read_csv(args.sim)
    df.sort_values(by="P", inplace=True)
    P = df["P"].to_numpy()
    v = df["v"].to_numpy()
    ax.plot(P, v, '-o')

    ax.set_xlabel(r"$P$[$\mu$m]")
    ax.set_ylabel(r"$v$[$\mu$m/s]")

    plt.show()


if __name__ == '__main__':
    main(sys.argv[1:])
