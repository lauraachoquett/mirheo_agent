#!/usr/bin/env python

import argparse
import numpy as np
import sdf_tools
import sdf_tools.Sdf as sdf
import sys

def create_segments(L: float=8,
                    R: float=1,
                    nx: int=3,
                    ny: int=3,
                    theta_deg: float=60, # in degrees
                    margin: float=2):

    assert R > 0
    assert L > 0
    assert theta_deg > 0
    assert theta_deg < 180
    assert nx > 1
    assert ny > 1

    theta = theta_deg * np.pi / 180
    t = theta/2

    dx = L * np.cos(t)
    dy = L * np.sin(t)

    Lx = (nx-1) * dx
    Ly = (ny-1) * dy
    Lz = 2 * (R + margin)

    xfrom = []
    yfrom = []
    xto = []
    yto = []

    for iy in range(ny):
        xfrom.append(0)
        yfrom.append(iy * dy)
        xto.append(min([Lx, iy*dy / np.tan(t)]))
        yto.append(0)

        xfrom.append(0)
        yfrom.append(iy * dy)
        xto.append(min([Lx, (Ly-iy*dy) / np.tan(t)]))
        yto.append(Ly)

    for ix in range(1, nx-1):
        xfrom.append(ix * dx)
        yfrom.append(0)
        xto.append(min([Lx, ix*dx + Ly / np.tan(t)]))
        yto.append(min([Ly, (Lx - ix*dx) * np.tan(t)]))

        xfrom.append(ix * dx)
        yfrom.append(Ly)
        xto.append(min([Lx, ix*dx + Ly / np.tan(t)]))
        yto.append(Ly - min([Ly, (Lx - ix*dx) * np.tan(t)]))


    nsegments = len(xfrom)
    x = np.zeros((nsegments,2))
    y = np.zeros((nsegments,2))
    z = np.zeros((nsegments,2)) + Lz/2

    x[:,0] = xfrom
    x[:,1] = xto

    y[:,0] = yfrom
    y[:,1] = yto

    r = np.full_like(xfrom, R)

    return x, y, z, r, (Lx, Ly, Lz)


def create_grid(L: float,
                R: float,
                nx: int,
                ny: int,
                theta_deg: float,
                margin: float,
                h: float,
                verbose: bool=False):

    x, y, z, r, domain = create_segments(L=L,
                                         R=R,
                                         nx=nx,
                                         ny=ny,
                                         theta_deg=theta_deg,
                                         margin=margin)

    nsegments = len(r)

    offs = [0, 0, 0]
    exts = domain
    dims = [max([1, int(d/h)]) for d in domain]

    if verbose:
        print(f"Grid dimensions: {dims}")
        print(f"Number of segments: {nsegments}")
        print(f"L, R = {L}, {R}")

    for i in range(nsegments):
        capsule = sdf.Capsule(start  = [x[i,0], y[i,0], z[i,0]],
                              end    = [x[i,1], y[i,1], z[i,1]],
                              radius = r[i],
                              inside = True)

        if i == 0:
            my_sdf = capsule
        else:
            my_sdf = sdf.SmoothUnion(my_sdf, capsule, k=R/24)


    grid = sdf_tools.Grid.Uniform(dims, offs, exts)
    grid.evaluate_sdf(my_sdf)

    return grid, domain




def main(argv: list):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=str, default="out.bov", help="Output file.")
    parser.add_argument('--R', type=float, default=1, help="Radius of the pipes.")
    parser.add_argument('--L', type=float, default=16, help="Length of the pipe between two bifurcations.")
    parser.add_argument('--nx', type=int, default=3, help="Number of bifurcation per x period.")
    parser.add_argument('--ny', type=int, default=3, help="Number of bifurcation per y period.")
    parser.add_argument('--theta', type=float, default=60, help="Angle between the pipes, in degrees.")
    parser.add_argument('--margin', type=float, default=2, help="Margin in the y and z directions.")
    parser.add_argument('--h', type=float, default=0.25, help="Resolution when creating the sdf on grid.")
    args = parser.parse_args(argv)

    grid, domain = create_grid(L=args.L,
                               R=args.R,
                               nx=args.nx,
                               ny=args.ny,
                               theta_deg=args.theta,
                               margin=args.margin,
                               h=args.h)

    grid.dump(args.out)


if __name__ == '__main__':
    main(sys.argv[1:])
