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


def main(basedir: str,
         obj: str,
         out: str):
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

    volume_domain = params.L * np.pi * params.R**2

    if obj == 'ABF':
        volume_obj = params.mesh_ABF.volume
    else:
        volume_obj = params.mesh_ini.volume

    print(f"volume fraction: {num_OBJs * volume_obj / volume_domain} ({num_OBJs} {obj}s)")

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

    if out is not None:
        data = {'t': t}
        for i in range(num_OBJs):
            data[f"x{i}"] = xall[i,:]
            data[f"y{i}"] = yall[i,:]
            data[f"z{i}"] = zall[i,:]
        df = pd.DataFrame(data)
        df.to_csv(out, index=False)
    else:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        nl = 50
        for c in np.linspace(0, 1, nl):
            ym = np.quantile(yall, q=0+c/2, axis=0)
            yp = np.quantile(yall, q=1-c/2, axis=0)
            ax.fill_between(t, ym, yp, alpha=1/nl, color='royalblue', lw=0)

        #ax.plot(t, y, '-k')
        ax.plot(t, np.median(yall, axis=0), '-k')
        ax.set_xlabel(r'$t$ [s]')
        ax.set_ylabel(r'$y/R$')
        ax.set_xlim(np.min(t), np.max(t))
        ax.set_ylim(-1,1)
        plt.show()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('basedir', type=str, help="path of the simulation folder.")
    parser.add_argument('obj', type=str, choices=['ABF', 'RBC'], help="The kind of objects to measure.")
    parser.add_argument('--out', type=str, default=None, help="output csv file.")
    args = parser.parse_args()

    main(args.basedir, args.obj, args.out)
