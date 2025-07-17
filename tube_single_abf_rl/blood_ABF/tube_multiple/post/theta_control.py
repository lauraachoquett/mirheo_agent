#!/usr/bin/env python

import numpy as np
import os
import pandas as pd
import sys

here = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(here, ".."))

from cylinder.parameters import (SimulationParams,
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


def extract_trajectories(basedir: str,
                         obj: str):
    params = load_parameters(os.path.join(basedir, 'parameters.pkl'))
    objstats_fname = os.path.join(basedir, 'stats', obj+'.csv')
    df = pd.read_csv(objstats_fname)

    ids = df['objId'].to_numpy()
    comx = df['comx'].to_numpy()
    comy = df['comy'].to_numpy()
    comz = df['comz'].to_numpy()
    time = df['time'].to_numpy()

    num_OBJs = len(np.unique(ids))
    num_steps = len(ids) // num_OBJs

    xall = np.zeros((num_OBJs, num_steps))
    yall = np.zeros((num_OBJs, num_steps))
    zall = np.zeros((num_OBJs, num_steps))
    t = np.zeros(num_steps)

    for i in range(num_OBJs):
        idx = np.argwhere(ids == i).flatten()
        xall[i,:] = unroll_periodic(comx[idx], params.L)
        yall[i,:] = comy[idx]
        zall[i,:] = comz[idx]
        t = time[idx]

    cy = 2*params.rc + params.R
    cz = 2*params.rc + params.R
    yall -= cy
    zall -= cz
    t *= params.time_scale_.to('s').magnitude
    #l = params.length_scale_.to('um').magnitude
    l = 1/params.R
    xall *= l
    yall *= l
    zall *= l

    return t, xall, yall, zall


def extract_theta(basedir):
    s = basedir.split('_')
    i = s.index('theta')
    return float(s[i+1])


def main(basedirs: list,
         obj: str,
         out: str):

    all_t = []
    all_x = []
    all_y = []
    all_z = []
    theta = []

    for basedir in basedirs:
        t, x, y, z = extract_trajectories(basedir, obj)
        all_t.append(t)
        all_x.append(y)
        all_y.append(y)
        all_z.append(y)
        theta.append(extract_theta(basedir))

    n = min([len(t) for t in all_t])

    all_x = np.array([x[:,:n] for x in all_x])
    all_y = np.array([y[:,:n] for y in all_y])
    all_z = np.array([z[:,:n] for z in all_z])

    theta = np.array(theta)
    idx = np.argsort(theta)
    theta = theta[idx]

    all_x = all_x[idx,:,:]
    all_y = all_y[idx,:,:]
    all_z = all_z[idx,:,:]

    t = all_t[0][:n]

    i = -1
    y_median = np.median(all_y[:,:,i], axis=1)

    print(f"t = {t[i]} s")

    if out is None:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        nl=50
        for c in np.linspace(0, 1, nl):
            qm = (1-c)/2
            qp = 1-qm
            y_m = np.quantile(all_y[:,:,i], q=qm, axis=1)
            y_p = np.quantile(all_y[:,:,i], q=qp, axis=1)
            ax.fill_between(theta, y_m, y_p, lw=0, alpha=1/nl, color='royalblue')
        ax.plot(theta, y_median, '-k')
        ax.set_xlabel(r'$\theta$')
        ax.set_ylabel(r'$y/R$')
        ax.set_xlim(0, 90)
        ax.set_ylim(-1,1)
        plt.show()
    else:
        num_thetas, num_objs, _ = all_y.shape
        data = {'theta': theta}

        for j in range(num_objs):
            data[f"x{j}"] = all_x[:,j,i]
            data[f"y{j}"] = all_y[:,j,i]
            data[f"z{j}"] = all_z[:,j,i]
        df = pd.DataFrame(data)
        df.to_csv(out, index=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('basedirs', type=str, nargs='+', help="path of the simulation folders.")
    parser.add_argument('obj', type=str, choices=['ABF', 'RBC'], help="The kind of objects to measure.")
    parser.add_argument('--out', type=str, default=None, help="output csv file.")
    args = parser.parse_args()

    main(args.basedirs, args.obj, args.out)
