#!/usr/bin/env python

import argparse
import numpy as np
import sdf_tools
import sdf_tools.Sdf as sdf
from svg.path import parse_path
from xml.dom.minidom import parse


def extract_segments_from_svg(svg_fname: str,
                              dx: float=1.0):

    xfrom = []
    yfrom = []
    xto = []
    yto = []

    doc = parse(svg_fname)
    svg = doc.getElementsByTagName('svg')[0]
    path_elements = doc.getElementsByTagName('path')
    for path_element in path_elements:
        pass_desc = path_element.getAttributeNode('d')
        path = parse_path(pass_desc.value)
        n = int(1+path.length()/dx)

        points = [path.point(x) for x in np.linspace(0, 1, n, endpoint=True)]

        for i in range(n-1):
            xfrom.append(points[i].real)
            yfrom.append(points[i].imag)
            xto.append(points[i+1].real)
            yto.append(points[i+1].imag)

    return xfrom, yfrom, xto, yto


def segments_to_sdf(xfrom: np.ndarray,
                    yfrom: np.ndarray,
                    xto: np.ndarray,
                    yto: np.ndarray,
                    r: np.ndarray,
                    margin: float,
                    h: float):

    nsegments = len(r)
    x = np.zeros((nsegments, 2))
    y = np.zeros((nsegments, 2))
    z = np.zeros((nsegments, 2))
    x[:,0] = xfrom
    x[:,1] = xto
    y[:,0] = yfrom
    y[:,1] = yto
    z[:,0] = np.zeros_like(r)
    z[:,1] = np.zeros_like(r)


    offs = [np.min(x),
            np.min(y - r[:,np.newaxis]) - margin,
            np.min(z - r[:,np.newaxis]) - margin]

    exts = [np.ptp(x),
            np.max(y + r[:,np.newaxis]) - np.min(y - r[:,np.newaxis]) + 2*margin,
            np.max(z + r[:,np.newaxis]) - np.min(z - r[:,np.newaxis]) + 2*margin]


    dims = [int(l/h) for l in exts]
    print(f"dimensions of the grid: {dims}")

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

    return x, y, z, offs, exts, grid


def create_grid(*,
                svg_fname: str,
                radius: float,
                margin: float,
                seg_dx: float=1.0,
                grid_dx: float=0.5):

    xfrom, yfrom, xto, yto = extract_segments_from_svg(svg_fname=svg_fname,
                                                       dx=seg_dx)
    r = np.full_like(xfrom, radius)

    x, y, z, offs, exts, grid = segments_to_sdf(xfrom=xfrom,
                                                yfrom=yfrom,
                                                xto=xto,
                                                yto=yto,
                                                r=r,
                                                margin=margin,
                                                h=grid_dx)

    x += offs[0]
    y += offs[1]
    z += offs[2]
    domain = exts

    return grid, domain


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('svg', type=str, help="Input svg file.")
    parser.add_argument('out', type=str, help="Output sdf file (bov, sdf or vtk).")
    parser.add_argument('--R', type=float, default=2, help="Radius.")
    parser.add_argument('--resolution', type=float, default=0.5, help="grid resolution.")
    parser.add_argument('--margin', type=float, default=2, help="Additional margin in the y and z direction.")
    args = parser.parse_args()

    grid, domain = create_grid(svg_fname=args.svg,
                               radius=args.R,
                               margin=args.margin,
                               seg_dx=args.R/2,
                               grid_dx=args.resolution)

    grid.dump(args.out)
