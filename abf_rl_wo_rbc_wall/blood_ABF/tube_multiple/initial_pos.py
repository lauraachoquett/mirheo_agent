#!/usr/bin/env python

import numpy as np
import sys




class InitialPositions:
    def __init__(self, arg: str):
        args = arg.split(',')
        if len(args) < 1:
            raise ValueError(f"need at least one argument")

        self.desc = args[0]

        if self.desc == 'line':
            self.n = int(args[1])
        elif self.desc == 'rnd':
            self.n = int(args[1])
            self.seed = int(args[2])
        elif self.desc == 'packed':
            self.d = float(args[1])
        else:
            raise ValueError(f"unknown description {self.desc}")


    def get_com_q(self,
                  *,
                  L: float,
                  R: float,
                  domain: tuple,
                  ABF_extents: np.ndarray):

        assert len(domain) == 3
        assert len(ABF_extents) == 3

        def get_packed(d: float,
                       rnd_xoffset: bool=False):

            ABF_R = np.min(ABF_extents) / 2
            ABF_L = np.max(ABF_extents)

            dx = ABF_L + d
            a = 2 * ABF_R + d
            h = np.sqrt(3/4) * a

            nx = int(L / dx)
            ny = int(2 * R / a)
            nz = int(2 * R / h)

            dx = L / nx

            cy, cz = domain[1]/2, domain[2]/2

            com = []

            for iy in range(ny):
                for iz in range(nz):
                    y = cy + (iy - ny//2) * a
                    z = cz + (iz - nz//2) * h
                    if iz % 2 == 1:
                        y += a/2

                    if np.sqrt((y-cy)**2 + (z-cz)**2) + ABF_R + d < R:
                        offset = np.random.uniform(0,0.999) if rnd_xoffset else 1/2
                        for ix in range(nx):
                            x = (ix + offset) * dx
                            com.append([x,y,z])

            com_q = [x + [1, 0, 0, 0] for x in com]
            return com_q


        if self.desc == 'line':
            Dabf = ABF_extents[2]
            Dpipe = 2 * R
            n = self.n
            assert n * Dabf < Dpipe

            h = Dpipe / n

            com = [[domain[0]/4, domain[1]/2, domain[2]/2 - R + h/2 + i * h] for i in range(n)]
            com_q = [x + [1, 0, 0, 0] for x in com]

        elif self.desc == 'packed':
            com_q = get_packed(self.d)

        elif self.desc == 'rnd':
            d = 0.05 * np.min(ABF_extents)
            np.random.seed(self.seed)

            com_q = np.array(get_packed(d, rnd_xoffset=True))

            n = self.n
            if len(com_q) < n:
                raise RuntimeError(f"too few sites, got {len(com_q)}, required {n}")

            idx = np.arange(len(com_q))
            np.random.shuffle(idx)
            idx = idx[:n]

            com_q = com_q[idx,:].tolist()


        else:
            raise NotImplementedError(f"Not implemented for desc {self.desc}")

        return com_q






def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Run RBCs flowing in a capillary pipe.')
    parser.add_argument('--ABF-ic', type=InitialPositions, default='line,1', help="Initial configuration of ABFs.")
    args = parser.parse_args(argv)

    L = 100
    R = 50
    rc = 1
    domain = [L, 2 * (R + 2*rc), 2 * (R + 2*rc)]
    extents = [20, 8, 8]

    print(args.ABF_ic.get_com_q(L=L, R=R, domain=domain, ABF_extents=extents))

if __name__ == '__main__':
    main(sys.argv[1:])
