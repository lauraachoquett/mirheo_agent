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
    parser.add_argument('basedirs', nargs='+', type=str, help="Simulation directory(ies).")
    parser.add_argument('--out', type=str, default='GUI', help="The name of the output file.")
    parser.add_argument('--plane', type=str, default='tx', help="Projection plane in space time.")
    parser.add_argument('--out-csv', type=str, default=None, help="If set, will dump the data to that csv file.")
    args = parser.parse_args(argv)

    fig, ax = plt.subplots(figsize=[6.8, 6.8])

    out_data_csv = dict()

    for basedir in args.basedirs:
        csv = os.path.join(basedir, 'obj_stats', 'ABF.csv')
        params = os.path.join(basedir, 'parameters.pkl')

        p = load_parameters(params)
        df = pd.read_csv(csv)

        obj_ids = df['objId'].to_numpy()
        nswimmers = len(np.unique(obj_ids))

        def get_varname(axis: str):
            if axis == 't':
                return 'time'
            else:
                return 'com'+axis


        def convert(data: np.ndarray,
                    axis: str):
            if axis == 'x':
                data = unroll_periodic(data, p.Lx)

            if axis in "xyz":
                data /= p.ABF_length
            else:
                data *= p.omega

            return data

        def make_label(axis: str):
            if axis in "xyz":
                return f"{axis} / L"
            return "t \omega"

        for i in range(nswimmers):
            axis0 = args.plane[0]
            axis1 = args.plane[1]
            key0 = get_varname(axis0)
            key1 = get_varname(axis1)

            idx = np.argwhere(obj_ids == i).flatten()

            x = df[key0].to_numpy()[idx]
            y = df[key1].to_numpy()[idx]

            x = convert(x, axis0)
            y = convert(y, axis1)

            ax.plot(x, y)

            out_data_csv[f"{axis0}{i}"] = x
            out_data_csv[f"{axis1}{i}"] = y

    ax.set_xlabel(f'${make_label(axis0)}$')
    ax.set_ylabel(f'${make_label(axis1)}$')

    if args.out_csv is not None:
        df = pd.DataFrame(out_data_csv)
        df.to_csv(args.out_csv, index=False)

    if args.out == 'GUI':
        plt.show()
    else:
        plt.savefig(args.out, transparent=True)

if __name__ == '__main__':
    main(sys.argv[1:])
