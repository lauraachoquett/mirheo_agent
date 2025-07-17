#!/usr/bin/env python

import argparse
from dataclasses import dataclass
import numpy as np
import pickle
import pint
import sys

import dpdprops

@dataclass
class Parameters:
    bead_radius: float
    nd: float
    rc: float
    mass: float
    dpd_oo: dpdprops.DPDParams
    domain: tuple
    omega: float
    magn_B: float
    magn_m: tuple
    magn_permeability: float
    length_scale_: pint.Quantity
    time_scale_: pint.Quantity
    mass_scale_: pint.Quantity
    magn_scale_: pint.Quantity


def create_parameters(ureg: pint.UnitRegistry,
                      L_: pint.Quantity,
                      bead_radius: float,
                      bead_radius_: pint.Quantity,
                      magn_B_: pint.Quantity,
                      magn_m_: pint.Quantity,
                      frequency_: pint.Quantity,
                      Re: float,
                      Ma: float,
                      verbose: bool=True):

    assert frequency_.check('[frequency]')

    # physical quantities
    eta_ = 1.0016 * ureg.mPa * ureg.s
    rho_ = 1e3 * ureg.kg / ureg.m**3
    magn_permeability_ = 1.256627e-6 * ureg.H / ureg.m # magn permeability of water

    length_scale_ = (bead_radius_ / bead_radius).to(ureg.um)

    omega_ = 2 * np.pi * frequency_
    U_ = bead_radius_ * omega_
    Re_ = float(U_ * bead_radius_ * rho_ / eta_)

    torque_ratio_ = float(magn_B_ * magn_m_ / (bead_radius_**3 * eta_ * omega_))


    # simulation parameters
    nd   = 10
    mass = 1
    rc   = 1

    omega  = 0.1       # sets the time scale
    rho    = nd * mass # sets the mass scale
    magn_B = 10000     # sets the magnetic scale


    time_scale_   = (omega / omega_).to(ureg.s)
    mass_scale_   = (rho_/rho * length_scale_**3).to(ureg.kg)
    magn_scale_   = (magn_B_ / magn_B).to(ureg.T)

    # setup DPD parameters

    U = bead_radius * omega
    dpd_oo = dpdprops.create_dpd_params_from_Re_Ma(Re=Re,
                                                   Ma=Ma,
                                                   U=U,
                                                   L=bead_radius,
                                                   nd=nd,
                                                   rc=rc,
                                                   mass=mass)

    # setup magnetic parameters

    torque_ratio = torque_ratio_
    magn_m = torque_ratio * omega * bead_radius**3 * dpd_oo.dynamic_viscosity() / magn_B
    magn_m = np.array([0, magn_m, 0])

    magn_permeability = float(magn_permeability_ * mass_scale_ / (length_scale_ * time_scale_**2 * magn_scale_**2))

    L = float(L_ / length_scale_)

    if verbose:
        print(f"Re = {Re_} (sim: {Re})")
        print(f"domain: {L_} x {L_} x {L_} (sim: {L} x {L} x {L})")
        print(f"omega: {omega_} (sim: {omega})")
        print(f"length scale: {length_scale_}")
        print(f"time   scale: {time_scale_}")
        print(f"mass   scale: {mass_scale_}")
        print(f"magn   scale: {magn_scale_}")
        print(f"dpd parameters: {dpd_oo}")
        print(f"magn_B: {magn_B_} (sim: {magn_B})")
        print(f"magn_m: (sim: {magn_m})")
        print(f"magn_perm: {magn_permeability_} (sim: {magn_permeability})")

    return Parameters(bead_radius=bead_radius,
                      nd=nd,
                      rc=rc,
                      mass=mass,
                      dpd_oo=dpd_oo,
                      domain=(L, L, L),
                      omega=omega,
                      magn_B=magn_B,
                      magn_m=magn_m,
                      magn_permeability=magn_permeability,
                      length_scale_=length_scale_,
                      time_scale_=time_scale_,
                      mass_scale_=mass_scale_,
                      magn_scale_=magn_scale_)


def dump_parameters(p: Parameters,
                    filename: str):
    import pickle
    with open(filename, 'wb') as f:
        pickle.dump(p, f, pickle.HIGHEST_PROTOCOL)


def load_parameters(filename: str):
    import pickle
    with open(filename, 'rb') as f:
        return pickle.load(f)


def get_quantity(ureg: pint.UnitRegistry,
                 s: str):
    # we avoid "/" and spaces so that directory names are not weird
    if '/' in s or ' ' in s:
        raise ValueError(f"'{s}' should not contain '/' or spaces, otherwise it will create directories with undesired names. Use '_per_' and '_' instead.")
    return ureg(s.replace('_', ' '))


def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--Re', type=float, default=0.1, help="The simulation Reynolds number.")
    parser.add_argument('--Ma', type=float, default=0.01, help="The simulation Mach number.")
    parser.add_argument('--frequency', type=str, default="10_Hz", help="Rotation frequency of the magnetic field.")
    parser.add_argument('--magn-m', type=str, default="1e-14_N_m_per_T", help="Magnitude of the magnetic moment of the swimmer.")
    parser.add_argument('--L', type=str, default="100_um", help="Size of the periodic domain.")
    parser.add_argument('--R', type=str, default="5_um", help="Radius of a single bead.")
    parser.add_argument('--R-sim', type=float, default=5, help="Radius of a single bead, in simulation units.")
    parser.add_argument('--out-params', type=str, default="parameters.pkl", help="Save all simulation parameters to this file.")
    args = parser.parse_args(argv)

    ureg = pint.UnitRegistry()

    L_ = get_quantity(ureg, args.L)

    # magnetic field
    frequency_ = get_quantity(ureg, args.frequency)

    magn_B_ = 3 * ureg.mT
    magn_m_ = get_quantity(ureg, args.magn_m)

    bead_radius   = args.R_sim
    bead_radius_  = get_quantity(ureg, args.R)

    p = create_parameters(ureg = ureg,
                          L_=L_,
                          Re=args.Re,
                          Ma=args.Ma,
                          bead_radius=bead_radius,
                          bead_radius_=bead_radius_,
                          magn_B_ = magn_B_,
                          magn_m_ = magn_m_,
                          frequency_ = frequency_)

    dump_parameters(p, args.out_params)


if __name__ == '__main__':
    main(sys.argv[1:])
