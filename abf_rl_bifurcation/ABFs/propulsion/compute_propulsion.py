#!/usr/bin/env python

import glob
import numpy as np
import os
import pandas as pd
import re
import sys

def get_velocities(dirname: str):
    rexf = '[-+]?\d*\.\d+|\d+'
    m = re.findall(f"U_({rexf})_({rexf})_({rexf})_W_({rexf})_({rexf})_({rexf})", dirname)[0]

    U = [float(m[0]),
         float(m[1]),
         float(m[2])]

    W = [float(m[3]),
         float(m[4]),
         float(m[5])]

    return U, W


def get_average_FT(dirname: str):
    fname = os.path.join(dirname, 'rigid.csv')
    df = pd.read_csv(fname)

    t = df['time'].to_numpy()

    fx = df['fx'].to_numpy()
    fy = df['fy'].to_numpy()
    fz = df['fz'].to_numpy()

    Tx = df['Tx'].to_numpy()
    Ty = df['Ty'].to_numpy()
    Tz = df['Tz'].to_numpy()

    tmax = np.max(t)
    idx = np.argwhere(t > tmax/5)

    f = [np.mean(fx[idx]),
         np.mean(fy[idx]),
         np.mean(fz[idx])]
    T = [np.mean(Tx[idx]),
         np.mean(Ty[idx]),
         np.mean(Tz[idx])]

    return f, T


def get_params(dirname: str):

    fname = os.path.join(dirname, "params.txt")

    params = dict()

    with open(fname, "r") as f:
        for line in f.readlines():
            rexf = '[-+]?\d*\.\d+|\d+'
            m = re.findall(f"([a-z]+) = ({rexf})", line)[0]
            params[m[0]] = float(m[1])

    return params


def process_pin_case(dirname: str, M: dict):
    U, W = get_velocities(dirname)
    f, T = get_average_FT(dirname)
    p = get_params(dirname)

    eta = p['eta']
    l = p['l']

    # check that only one component is non zero
    assert np.max(np.abs(U + W)) == np.sum(np.abs(U + W))

    i = np.argmax(np.abs(U + W))

    if i < 3:
        Ai = -U[i] / f[i]
        # make it dimless
        Ai = Ai * eta * l
        print(f"A{i} = {Ai}")
        M[f'A{i}'] = Ai
        # if i == 0:
        #     Bi = U[i] / T[i]
        #     Bi = Bi * eta * l**2
        #     print(f"B{i} = {Bi}")

    else:
        # if i == 3:
        #     Bi = W[i-3] / f[i-3]
        #     Bi = Bi * eta * l**2
        #     print(f"B{i-3} = {Bi}")

        Ci = -W[i-3] / T[i-3]
        Ci = Ci * eta * l**3
        print(f"C{i-3} = {Ci}")
        M[f'C{i-3}'] = Ci


def compute_mean_velocity(dirname: str):
    fname = os.path.join(dirname, 'rigid.csv')
    df = pd.read_csv(fname)
    t = df['time'].to_numpy()
    x = df['comx'].to_numpy()
    p = np.polyfit(t, x, deg=1)
    return p[0]


def process_forward_case(dirname: str, M: dict):
    U, W = get_velocities(dirname)
    assert W[0] > 0
    assert U[0] < 0

    U = compute_mean_velocity(dirname)
    p = get_params(dirname)
    l = p['l']

    U = U/l
    W = W[0]

    B0 = U / W * M["C0"]
    print(f"B0 = {B0}")
    M['B0'] = B0




def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('base_dir', type=str, help="Base directory that contains the simulations output.")
    args = parser.parse_args(argv)

    M = dict()

    cases = sorted(glob.glob(os.path.join(args.base_dir, "pin_U_*_*_*_W_*_*_*")))[::-1]

    for c in cases:
        process_pin_case(c, M)

    cases = glob.glob(os.path.join(args.base_dir, "stats_U_*_*_*_W_*_*_*"))

    for c in cases:
        process_forward_case(c, M)

    B = 1 # 1mT
    Vmax = 1 # body lengths per second
    m = Vmax / M['B0'] / B
    print(f"m = {m}")

    wc = M['C0'] * m * B
    print(f"wc = {wc}")

if __name__ == '__main__':
    main(sys.argv[1:])
