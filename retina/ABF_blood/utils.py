import numpy as np
import os
import pyvista as pv
from scipy.interpolate import interp1d

def coordinate_in_path_ref_3D(p, x, t, n, b):
    B = np.column_stack((t, n, b))
    x = np.array(x) - np.array(p)
    return B.T @ x


def coordinate_in_global_ref_3D(p, x, t, n, b):
    B = np.column_stack((t, n, b))
    x = B @ x + p
    return x

def generate_simple_line(p_0, p_target, nb_points):
    t = np.linspace(0, 1, nb_points)
    path = p_0 * t[:, None] + (1 - t)[:, None] * p_target
    d = np.linalg.norm((p_target - p_0))
    return np.flip(path, axis=0), d

def generate_helix(num_points=1000, radius=1.0, pitch=0.1, turns=3,
                   clockwise=True, x_0=0.0, y_0=0.0, z_0=0.0):
    t = np.linspace(0, 2 * np.pi * turns, num_points)

    if clockwise:
        z = radius * np.sin(t)
        y = radius * np.cos(t)
    else:
        z = radius * np.cos(t)
        y = radius * np.sin(t)

    dz = z_0 - z[0]
    dy = y_0 - y[0]
    dx = x_0 - pitch * t[0] / (2 * np.pi)

    z += dz
    y += dy
    x = pitch * t / (2 * np.pi) + dx

    curve = np.stack((x, y, z), axis=1)
    return curve

def paraview_export(path_physical,  output_save_path):
    os.makedirs(output_save_path,exist_ok=True)
    if path_physical is not None and len(path_physical) > 0:
        print(f"Path points shape: {path_physical.shape}")
        if path_physical.shape[1] == 3:
            path_polydata = pv.PolyData(path_physical)
            
            if len(path_physical) > 1:
                lines = []
                for i in range(len(path_physical) - 1):
                    lines.extend([2, i, i + 1])  
                path_polydata.lines = np.array(lines)
                print(f"Nombre de segments créés: {len(path_physical) - 1}")
            
            path_output = os.path.join(output_save_path,'path.vtp')
            path_polydata.save(path_output)
            print(f"Chemin sauvegardé: {path_output}")
            

import numpy as np

def normalize(v):
    return v / np.linalg.norm(v)

def double_reflection_rmf(points):
    points = np.array(points)
    n = len(points)

    # Tangents
    tangents = [normalize(points[i+1] - points[i]) for i in range(n - 1)]
    
    # Initial frame: choose arbitrary normal N0 orthogonal to T0
    T0 = tangents[0]
    arbitrary = np.array([0, 0, 1])
    if np.allclose(np.cross(T0, arbitrary), 0):
        arbitrary = np.array([1, 0, 0])
    N0 = normalize(np.cross(np.cross(T0, arbitrary), T0))
    B0 = np.cross(T0, N0)

    frames = [(T0, N0, B0)]
    N_prev = N0

    for i in range(1, n - 1):
        v1 = tangents[i - 1]
        v2 = tangents[i]

        # First reflection
        r = v2 + v1
        if np.linalg.norm(r) < 1e-10:  # 180° turn: reset
            N_curr = N_prev
        else:
            r = normalize(r)
            N_prime = N_prev - 2 * np.dot(N_prev, r) * r

            # Second reflection
            r = v2
            N_curr = N_prime - 2 * np.dot(N_prime, r) * r

        B_curr = np.cross(v2, N_curr)
        frames.append((v2, N_curr, B_curr))
        N_prev = N_curr
    frames.append((v2, N_curr, B_curr))
    # Convert frames to arrays
    T = np.array([f[0] for f in frames])
    N = np.array([f[1] for f in frames])
    B = np.array([f[2] for f in frames])
    return T, N, B

def resample_path(path, n_points=500):
    distances = np.sqrt(np.sum(np.diff(path, axis=0) ** 2, axis=1))
    cumulative = np.concatenate(([0], np.cumsum(distances)))
    total_length = cumulative[-1]
    fx = interp1d(cumulative, path[:, 0], kind="linear")
    fy = interp1d(cumulative, path[:, 1], kind="linear")
    fz = interp1d(cumulative, path[:, 2], kind="linear")
    new_distances = np.linspace(0, total_length, n_points)
    return np.stack((fx(new_distances), fy(new_distances),fz(new_distances)), axis=1), np.sum(distances)



def sdf_read(path):
    with open(path, "rb") as f:
        L = [float(l) for l in f.readline().decode().split()]
        N = [int(n) for n in f.readline().decode().split()]
        field = np.fromfile(f, dtype=np.float32).reshape(N[::-1])
    h = [l/n for l, n in zip(L, N)]
    return N, h, field


def sdf_write(path, field, shape, space):
    L = np.array(space) * np.array(shape)

    with open(path, "w") as f:
        print("{} {} {}".format(*L), file=f)
        print("{} {} {}".format(*shape), file=f)
        data = np.float32(field)
        data.tofile(f)


def force_read(path):
    with open(path, "rb") as f:
        L = [float(l) for l in f.readline().decode().split()]
        N = [int(n) for n in f.readline().decode().split()]
        field = np.fromfile(f, dtype=np.float32).reshape((*N[::-1], 4))
    h = [l/n for l, n in zip(L, N)]
    return N, h, field


def force_write(path, field, shape, space):
    L = np.array(space) * np.array(shape)

    with open(path, "w") as f:
        print("{} {} {}".format(*L), file=f)
        print("{} {} {}".format(*shape), file=f)
        data = np.float32(field)
        data.tofile(f)



if __name__ == '__main__':
    R = 17.96
    pos = []
    path = generate_helix(num_points=1000, radius=R*3/4, pitch=3, turns=4, clockwise=True,x_0 = pos[0],y_0=pos[1],z_0=pos[2])
    paraview_export(path,  'paraview_export/')
