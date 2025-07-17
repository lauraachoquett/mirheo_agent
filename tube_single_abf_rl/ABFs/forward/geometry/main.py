#!/usr/bin/env python

import numpy as np
import trimesh

def create_mesh(L: float=15,
                D: float=5,
                thickness: float=1,
                pitch: float=3,
                dx: float=0.2):

    nturns = (L - 2*thickness) / pitch
    T = nturns * 2 * np.pi
    a = (D - 2*thickness) / 2
    b = (L - 2*thickness) / T

    ax = thickness/2
    ay = thickness/2

    length = T * np.sqrt(a**2 + b**2)
    nt = int(length / dx)
    ntheta = int(np.pi * thickness / dx)

    if False:
        print(f"length: {length}")
        print(f"curvature: {abs(a) / (a**2 + b**2)}")
        print(f"torsion: {b / (a**2 + b**2)}")

    vertices = []
    ts     = np.linspace(0, T, nt, endpoint=True)
    thetas = np.linspace(0, 2 * np.pi, ntheta, endpoint=True)

    for t in ts:
        r0 = np.array([a * np.cos(t),
                       a * np.sin(t),
                       b * t])

        N = np.array([-np.cos(t),
                      -np.sin(t),
                      0])

        B = np.array([+b * np.sin(t),
                      -b * np.cos(t),
                      +a]) / np.sqrt(a**2 + b**2)

        s = t / T
        delta = 0.02
        f = 1

        if s < delta:
            x = s / delta
            f = np.sqrt(1 - (1-x)**2)
        elif s > 1 - delta:
            x = (1-s) / delta
            f = np.sqrt(1 - (1-x)**2)

        ax_ = f * ax
        ay_ = f * ay

        for theta in thetas:
            dr = ax_ * np.cos(theta) * N + ay_ * np.sin(theta) * B
            r = r0 + dr
            vertices.append(r)


    k = 0
    faces = []
    for i in range(nt-1):
        for j in range(ntheta-1):
            faces.append((k + j, k + j + 1, k + ntheta + j))
            faces.append((k + j + 1, k + ntheta + j + 1, k + ntheta + j))
        k += ntheta

    mesh = trimesh.Trimesh(vertices=vertices,
                           faces=faces,
                           process=True) # to merge vertices

    faces_cleaned = []
    for i, f in enumerate(mesh.faces):
        if not (f[0] == f[1] or f[1] == f[2] or f[0] == f[2]):
            faces_cleaned.append(f)

    mesh = trimesh.Trimesh(vertices=mesh.vertices,
                           faces=faces_cleaned,
                           process=True)

    assert mesh.is_watertight

    return mesh


if __name__ == '__main__':

    for p in range(2, 10):
        mesh = create_mesh(L=15, D=5, thickness=1, pitch=p)
        mesh.export(f"ABF_P_{p}.ply")
