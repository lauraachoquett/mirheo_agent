#!/usr/bin/env python

import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import re
import sys

def quat_conjugate(q: np.ndarray):
    assert len(q.shape) == 2
    assert q.shape[1] == 4
    res = np.zeros_like(q)
    res[:,0] = q[:,0]
    res[:,1:]  = -q[:,1:]
    return res

def quat_product(a: np.ndarray, b: np.ndarray):
    assert a.shape == b.shape
    assert len(a.shape) == 2
    assert a.shape[1] == 4

    aw = a[:,0]; ax = a[:,1]; ay = a[:,2]; az = a[:,3]
    bw = b[:,0]; bx = b[:,1]; by = b[:,2]; bz = b[:,3]

    res = np.zeros_like(a)
    res[:,0] = aw * bw - ax * bx - ay * by - az * bz
    res[:,1] = aw * bx + ax * bw + ay * bz - az * by
    res[:,2] = aw * by - ax * bz + ay * bw + az * bx
    res[:,3] = aw * bz + ax * by - ay * bx + az * bw

    return res

def compute_wobbling_angles(obj_stats_csv: str):

    df = pd.read_csv(obj_stats_csv)

    obj_ids = df['objId'].to_numpy()
    nswimmers = len(np.unique(obj_ids))
    assert nswimmers == 1

    qw = df['qw'].to_numpy()
    qx = df['qx'].to_numpy()
    qy = df['qy'].to_numpy()
    qz = df['qz'].to_numpy()

    q = np.concatenate((qw.reshape((-1,1)),
                        qx.reshape((-1,1)),
                        qy.reshape((-1,1)),
                        qz.reshape((-1,1))),
                       axis=1)

    ex = np.zeros_like(q)
    ex[:,1] = 1

    qstar = quat_conjugate(q)

    n = quat_product(quat_product(qstar, ex), q)
    # normalize
    n /= np.sqrt(np.sum((n*n)[:,1:], axis=1))[:,np.newaxis]

    costhetas = np.sum((ex * n)[:,1:], axis=1)
    thetas = np.arccos(costhetas)
    return thetas

def get_freq(basedir: str):
    rexf = '[-+]?\d*\.\d+|\d+'
    matches = re.findall(f"freq_({rexf})_Hz", basedir)
    assert len(matches) == 1
    return float(matches[0])



def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('paths', type=str, nargs='+', help="path of the simulation folder(s).")
    parser.add_argument('--out', type=str, default="results/f_theta.csv", help="output csv file.")
    args = parser.parse_args(argv)

    all_f = list()
    all_thetas = list()

    for basedir in args.paths:
        obj_stats_csv = os.path.join(basedir, 'obj_stats', 'ABF.csv')

        thetas = compute_wobbling_angles(obj_stats_csv)

        if 1:
            print(180*thetas/np.pi)
        if 0:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            ax.plot(range(len(thetas)), thetas * 180/np.pi)
            ax.set_ylabel(r"$\theta$")
            ax.set_xlabel("time stamp")
            plt.show()

        theta = np.mean(thetas[:len(thetas)//2])
        f = get_freq(basedir)
        all_f.append(f)
        all_thetas.append(theta)
        print(f"f = {f} Hz, theta = {theta}")

    df = pd.DataFrame({"f" : all_f,
                       "theta": all_thetas})
    df.to_csv(args.out, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
