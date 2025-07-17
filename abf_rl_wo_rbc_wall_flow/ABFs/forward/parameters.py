#!/usr/bin/env python

import argparse
from dataclasses import dataclass
import numpy as np
import pickle
import pint
import trimesh
import sys

from geometry import create_mesh
import dpdprops

@dataclass
class Parameters:
    ABF_mesh: trimesh.Trimesh
    ABF_length: float
    ABF_radius: float
    ABF_thickness: float
    ABF_pitch: float
    nd: float
    rc: float
    mass: float
    dpd_oo: dpdprops.DPDParams
    domain: tuple
    omega: float
    magn_B: float
    magn_m: tuple
    length_scale_: pint.Quantity
    time_scale_: pint.Quantity
    mass_scale_: pint.Quantity
    magn_scale_: pint.Quantity


def create_parameters(ureg: pint.UnitRegistry,
                      Lx_: pint.Quantity,
                      Ly_: pint.Quantity,
                      Lz_: pint.Quantity,
                      ABF_length: float,
                      ABF_length_: pint.Quantity,
                      ABF_pitch_: pint.Quantity,
                      ABF_radius_: pint.Quantity,
                      ABF_thickness_: pint.Quantity,
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

    # ABF geometry
    length_scale_ = (ABF_length_ / ABF_length).to(ureg.um)
    ABF_pitch = float(ABF_pitch_ / length_scale_)
    ABF_radius = float(ABF_radius_ / length_scale_)
    ABF_thickness = float(ABF_thickness_ / length_scale_)

    omega_ = 2 * np.pi * frequency_
    U_ = ABF_radius_ * omega_
    Re_ = float(U_ * ABF_radius_ * rho_ / eta_)

    torque_ratio_ = float(magn_B_ * magn_m_ / (ABF_radius_**3 * eta_ * omega_))


    # simulation parameters
    nd   = 10
    mass = 1
    rc   = 1

    omega  = 0.1       # sets the time scale
    rho    = nd * mass # sets the mass scale
    magn_B = 10        # sets the magnetic scale


    time_scale_   = (omega / omega_).to(ureg.s)
    mass_scale_   = (rho_/rho * length_scale_**3).to(ureg.kg)
    magn_scale_   = (magn_B_ / magn_B).to(ureg.T)


    ABF_mesh = create_mesh(L=ABF_length,
                           D=2*ABF_radius,
                           thickness=ABF_thickness,
                           pitch=ABF_pitch,
                           dx=ABF_thickness * np.pi/16)

    # align the mesh to its principal axes
    I = ABF_mesh.moment_inertia

    Lambda, W = trimesh.inertia.principal_axis(I)
    vertices = np.array(ABF_mesh.vertices)
    vertices -= ABF_mesh.center_mass
    vertices = np.dot(vertices, W.T)
    ABF_mesh.vertices = vertices
    if ABF_mesh.volume < 0:
        ABF_mesh.faces = ABF_mesh.faces[:,::-1]

    # setup DPD parameters

    U = ABF_radius * omega
    dpd_oo = dpdprops.create_dpd_params_from_Re_Ma(Re=Re,
                                                  Ma=Ma,
                                                  U=U,
                                                  L=ABF_radius,
                                                  nd=nd,
                                                  rc=rc,
                                                  mass=mass)

    # setup magnetic parameters

    torque_ratio = torque_ratio_
    magn_m = torque_ratio * omega * ABF_radius**3 * dpd_oo.dynamic_viscosity() / magn_B
    magn_m = np.array([0, magn_m, 0])
    magn_m = np.dot(magn_m, W.T)
    if magn_m[2] < 0:
        magn_m *= -1
    magn_m = tuple(magn_m)


    Lx = float(Lx_ / length_scale_)
    Ly = float(Ly_ / length_scale_)
    Lz = float(Lz_ / length_scale_)

    if verbose:
        print(f"Re = {Re_} (sim: {Re})")
        print(f"domain: {Lx_} x {Ly_} x {Lz_} (sim: {Lx} x {Ly} x {Lz})")
        print(f"omega: {omega_} (sim: {omega})")
        print(f"length scale: {length_scale_}")
        print(f"time   scale: {time_scale_}")
        print(f"mass   scale: {mass_scale_}")
        print(f"magn   scale: {magn_scale_}")
        print(f"dpd parameters: {dpd_oo}")
        print(f"magn_B: {magn_B_} (sim: {magn_B})")
        print(f"magn_m: (sim: {magn_m})")

    return Parameters(ABF_mesh=ABF_mesh,
                      ABF_length=ABF_length,
                      ABF_radius=ABF_radius,
                      ABF_pitch=ABF_pitch,
                      ABF_thickness=ABF_thickness,
                      nd=nd,
                      rc=rc,
                      mass=mass,
                      dpd_oo=dpd_oo,
                      domain=(Lx, Ly, Lz),
                      omega=omega,
                      magn_B=magn_B,
                      magn_m=magn_m,
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
    parser.add_argument('--L', type=str, nargs=3, default=["200_um", "100_um", "100_um"], help="Size of the periodic domain.")
    parser.add_argument('--ABF-length', type=str, default="15_um", help="Length of the helix.")
    parser.add_argument('--ABF-pitch', type=str, default="4_um", help="Pitch of the helix.")
    parser.add_argument('--ABF-radius', type=str, default="10_um", help="Radius of the helix.")
    parser.add_argument('--ABF-thickness', type=str, default="10_um", help="Thickness of the helix.")
    parser.add_argument('--ABF-L-sim', type=float, default=20, help="Length of the ABF, in simulation units.")
    parser.add_argument('--out-params', type=str, default="parameters.pkl", help="Save all simulation parameters to this file.")
    parser.add_argument('--out-mesh', type=str, default="ABF_mesh.ply", help="Save the ABF mesh to this file.")

    args = parser.parse_args(argv)

    ureg = pint.UnitRegistry()

    Lx_, Ly_, Lz_ = [get_quantity(ureg, L) for L in args.L]

    # magnetic field
    frequency_ = get_quantity(ureg, args.frequency)

    magn_B_ = 3 * ureg.mT
    magn_m_ = get_quantity(ureg, args.magn_m)

    ABF_length     = args.ABF_L_sim
    ABF_length_    = get_quantity(ureg, args.ABF_length)
    ABF_pitch_     = get_quantity(ureg, args.ABF_pitch)
    ABF_radius_    = get_quantity(ureg, args.ABF_radius)
    ABF_thickness_ = get_quantity(ureg, args.ABF_thickness)

    p = create_parameters(ureg = ureg,
                          ABF_length=ABF_length,
                          Lx_=Lx_, Ly_=Ly_, Lz_=Lz_,
                          Re=args.Re,
                          Ma=args.Ma,
                          ABF_length_  = ABF_length_,
                          ABF_pitch_  = ABF_pitch_,
                          ABF_radius_ = ABF_radius_,
                          ABF_thickness_ = ABF_thickness_,
                          magn_B_ = magn_B_,
                          magn_m_ = magn_m_,
                          frequency_ = frequency_)

    p.ABF_mesh.export(args.out_mesh)
    dump_parameters(p, args.out_params)


if __name__ == '__main__':
    main(sys.argv[1:])
