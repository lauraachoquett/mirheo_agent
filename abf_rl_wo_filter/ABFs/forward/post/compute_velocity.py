#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import re
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


def get_frequency(basedir: str):
    rexf = '[-+]?\d*\.\d+|\d+'
    matches = re.findall(f"freq_({rexf})_Hz", basedir)
    assert len(matches) == 1
    return float(matches[0])

def get_L(basedir: str):
    rexf = '[-+]?\d*\.\d+|\d+'
    matches = re.findall(f"L_({rexf})_um_({rexf})_um_({rexf})_um", basedir)
    assert len(matches) == 1
    return [float(m) for m in matches[0]]

def get_nondim_param(basedir: str, name: str):
    rexf = '[-+]?\d*\.\d+|\d+'
    matches = re.findall(f"{name}_({rexf})", basedir)
    assert len(matches) == 1
    return float(matches[0])

def get_varying_param_val(varname: str, basedir: str):
    if varname == 'freq':
        return get_frequency(basedir)
    elif varname == 'L':
        return get_L(basedir)[0]
    else:
        return get_nondim_param(basedir, varname)

def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('basedirs', type=str, nargs='+', help="Simulation directory.")
    parser.add_argument('--out', type=str, default=None, help="The name of the output csv file.")
    parser.add_argument('--against', type=str, default="freq", help="varying parameter.")
    args = parser.parse_args(argv)

    varying_param_vals = []
    velocities = []

    for basedir in args.basedirs:
        varying_param_val = get_varying_param_val(args.against, basedir)
        csv    = os.path.join(basedir, 'obj_stats', 'ABF.csv')
        params = os.path.join(basedir, 'parameters.pkl')

        p  = load_parameters(params)
        df = pd.read_csv(csv)

        t = df['time'].to_numpy()
        x = df['comx'].to_numpy()

        start = len(x) // 2
        end = -1

        if True:
            pol = np.polyfit(t[start:end], x[start:end], deg=1)
            v = abs(pol[0])
        else:
            v = abs((x[end] - x[start]) / (t[end] - t[start]))

        v = v * p.length_scale_ / p.time_scale_
        print(f"{args.against} = {varying_param_val}, v = {v}")

        varying_param_vals.append(varying_param_val)
        velocities.append(v.to('um/s').magnitude)

    if args.out is not None:
        df = pd.DataFrame({args.against: varying_param_vals, "v": velocities})
        df.sort_values(by=args.against)
        df.to_csv(args.out, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
