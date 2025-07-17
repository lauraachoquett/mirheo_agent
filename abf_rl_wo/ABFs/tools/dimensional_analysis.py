#!/usr/bin/env python
# Copyright 2020 ETH Zurich. All Rights Reserved.

"""
Find what non dimensional quantities are relevant.
"""

from pint import UnitRegistry
ureg = UnitRegistry()

torque = 1 * ureg.m * ureg.newton
force = 1 * ureg.newton
dyn_visc = 1 * ureg.Pa * ureg.s

# propulsion matrix components, without the viscosity scaling
#
#
# V = 1/eta * (A * F + B * T)
# W = 1/eta * (B * F + C * T)

A = 1 / ureg.m
B = 1 / ureg.m**2
C = 1 / ureg.m**3

V = (A * force + B * torque) / dyn_visc
W = (B * force + C * torque) / dyn_visc

assert V.check('[length]/[time]')
assert W.check('1/[time]')

L = 1 * ureg.m

# sound speed in the fluid
Cs = 1 * ureg.m / ureg.s
# mass density of the fluid
rho = 1 * ureg.kg / ureg.m**3

pt = ureg.pi_theorem({'L': L,
                      'eta': dyn_visc,
                      'Omega': W,
                      'Cs': Cs,
                      'rho': rho,
                      'A': A,
                      'B': B,
                      'C': C})
for p in pt:
    print(p)
