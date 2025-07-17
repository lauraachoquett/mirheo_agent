#!/usr/bin/env python

import argparse
import numpy as np
import sdf_tools
import sdf_tools.Sdf as sdf
import sys

def create_segments(Lmain: float=8,
                    Rmain: float=3,
                    Lbranch: float=16,
                    Rbranch: float=1.5,
                    theta: float=60, # in degrees
                    margin: float=2,
                    ds: float=0.25,
                    tc: float=0.5):

    assert Rmain > 0
    assert Rbranch > 0
    assert Rmain >= Rbranch
    assert theta > 0
    assert theta < 180

    theta *= np.pi / 180

    Lx = 2 * (Lmain + Lbranch)
    Ly = 2 * (Rbranch + np.tan(theta/2) * Lbranch + 2*margin)
    Lz = 2 * (Rmain + margin)

    nmain = int(Lmain / ds)
    nbranch = int(Lbranch / ds)
    ntot = 2 * nmain + 4 * nbranch

    x = np.zeros((ntot,2))
    y = np.zeros((ntot,2))
    z = np.zeros((ntot,2))
    r = np.zeros(ntot)


    z[:,:] = Lz/2

    t = np.linspace(0, 1, nmain+1, endpoint=True)

    x[:nmain,0] = t[:-1] * Lmain
    x[:nmain,1] = t[+1:] * Lmain
    y[:nmain,:] = Ly/2
    r[:nmain] = Rmain


    t = np.linspace(0, 1, nbranch+1, endpoint=True)

    x[nmain:nmain+nbranch,0] = Lmain + t[:-1] * Lbranch
    x[nmain:nmain+nbranch,1] = Lmain + t[+1:] * Lbranch

    y[nmain:nmain+nbranch,0] = Ly/2 + np.sin(np.pi/2 * t[:-1]) * np.tan(theta/2) * Lbranch
    y[nmain:nmain+nbranch,1] = Ly/2 + np.sin(np.pi/2 * t[+1:]) * np.tan(theta/2) * Lbranch

    r[nmain:nmain+nbranch] = Rbranch + (Rmain-Rbranch) * np.exp(-t[:-1]/tc)

    x[nmain+nbranch:nmain+2*nbranch,0] = Lmain + t[:-1] * Lbranch
    x[nmain+nbranch:nmain+2*nbranch,1] = Lmain + t[+1:] * Lbranch

    y[nmain+nbranch:nmain+2*nbranch,0] = Ly/2 - np.sin(np.pi/2 * t[:-1]) * np.tan(theta/2) * Lbranch
    y[nmain+nbranch:nmain+2*nbranch,1] = Ly/2 - np.sin(np.pi/2 * t[+1:]) * np.tan(theta/2) * Lbranch

    r[nmain+nbranch:nmain+2*nbranch] = Rbranch + (Rmain-Rbranch) * np.exp(-t[:-1]/tc)

    # symmetry
    x[ntot//2:,:] = Lx - x[:ntot//2,:]
    y[ntot//2:,:] = y[:ntot//2,:]
    r[ntot//2:] = r[:ntot//2]

    return x, y, z, r, (Lx, Ly, Lz)


def create_grid(Lmain: float,
                Rmain: float,
                Rbranch: float,
                Lbranch: float,
                theta: float,
                margin: float,
                tc: float,
                h: float,
                ds: float,
                verbose: bool=False):

    x, y, z, r, domain = create_segments(Lmain=Lmain,
                                         Rmain=Rmain,
                                         Lbranch=Lbranch,
                                         Rbranch=Rbranch,
                                         theta=theta,
                                         margin=margin,
                                         ds=ds,
                                         tc=tc)

    nsegments = len(r)

    offs = [0, 0, 0]
    exts = domain
    dims = [max([1, int(d/h)]) for d in domain]

    if verbose:
        print(f"Grid dimensions: {dims}")
        print(f"Number of segments: {nsegments}")

    for i in range(nsegments):
        capsule = sdf.Capsule(start  = [x[i,0], y[i,0], z[i,0]],
                              end    = [x[i,1], y[i,1], z[i,1]],
                              radius = r[i],
                              inside = True)

        if i == 0:
            my_sdf = capsule
        else:
            my_sdf = sdf.Union(my_sdf, capsule)


    grid = sdf_tools.Grid.Uniform(dims, offs, exts)
    grid.evaluate_sdf(my_sdf)

    return grid, domain


def main(argv: list):
    import argparse
    parser = argparse.ArgumentParser(description='Create the sdf representation of a bifurcation on a grid.')
    parser.add_argument('--out', type=str, default="bifurcation.sdf", help="Output file.")
    parser.add_argument('--Rmain', type=float, default=3, help="Radius of the main branch.")
    parser.add_argument('--Rbranch', type=float, default=1.5, help="Radius of the child branches.")
    parser.add_argument('--Lmain', type=float, default=12, help="Length of the main branch.")
    parser.add_argument('--Lbranch', type=float, default=12, help="Length of the child branches.")
    parser.add_argument('--theta', type=float, default=60, help="Angle between the child branches, in degrees.")
    parser.add_argument('--margin', type=float, default=2, help="Margin in the y and z directions.")
    parser.add_argument('--ds', type=float, default=0.25, help="Resolution when creating the segments.")
    parser.add_argument('--h', type=float, default=0.25, help="Resolution when creating the sdf on grid.")
    parser.add_argument('--tc', type=float, default=0.5, help="Decay length for the transitions between the main and the child branches.")
    args = parser.parse_args(argv)

    grid, domain = create_grid(Lmain=args.Lmain,
                               Rmain=args.Rmain,
                               Lbranch=args.Lbranch,
                               Rbranch=args.Rbranch,
                               theta=args.theta,
                               margin=args.margin,
                               ds=args.ds,
                               tc=args.tc,
                               h=args.h,
                               verbose=True)
    grid.dump(args.out)


if __name__ == '__main__':
    main(sys.argv[1:])
