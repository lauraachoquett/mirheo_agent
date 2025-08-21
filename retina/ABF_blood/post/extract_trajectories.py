#!/usr/bin/env python

import argparse
import numpy as np
import os
import pandas as pd

from parameters import (ContactParams,
                        Parameters,
                        load_parameters)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=str, help='path to simulation directory')
    parser.add_argument('out', type=str, help='path to the output csv file')
    parser.add_argument('--Lx', type=float, default=1.0, help='rescale to this domain length along x')
    args = parser.parse_args()

    sim_path = args.path
    p = load_parameters(os.path.join(sim_path, "parameters.pkl"))

    df = pd.read_csv(os.path.join(sim_path, "obj_stats", "abf.csv"))

    x = df["comx"].to_numpy()
    y = df["comy"].to_numpy()
    z = df["comz"].to_numpy()

    scale = args.Lx / p.domain[0]

    x *= scale
    y *= scale
    z *= scale

    data = {
        'x': x,
        'y': y,
        'z': z
    }

    df = pd.DataFrame(data)
    df.to_csv(args.out, index=False)
