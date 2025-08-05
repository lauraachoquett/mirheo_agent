#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

def plot(fname: str, skip: int=5):
    df = pd.read_csv(fname)

    fig, ax = plt.subplots()

    t = df['time'].to_numpy()[skip:]
    x = df['comx'].to_numpy()[skip:]

    t -= t[0]
    x -= x[0]

    p = np.polyfit(t, x, deg=1)
    print(p)

    ax.plot(t, x, 'o')
    ax.plot(t, np.polyval(p, t), '--')

    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$x$')

    plt.show()

def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('stats', type=str, help="stats file.")
    args = parser.parse_args(argv)

    plot(args.stats)

if __name__ == '__main__':
    main(sys.argv[1:])
