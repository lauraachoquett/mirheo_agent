import numpy as np


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
