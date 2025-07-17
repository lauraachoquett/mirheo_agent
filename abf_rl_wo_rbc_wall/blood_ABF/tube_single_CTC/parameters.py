#!/usr/bin/env python

from copy import deepcopy
from dataclasses import dataclass
import numpy as np
import pickle
import pint
import trimesh
import sys

import dpdprops

@dataclass
class CapilarryFlowParams:
    rc: float
    m: float
    nd: float
    RA: float
    kBT: float
    L: float
    R: float
    mesh_ini: trimesh.Trimesh
    mesh_ref: trimesh.Trimesh
    mesh_ABF: trimesh.Trimesh
    mesh_CTC: trimesh.Trimesh
    rbc_params: dpdprops.MembraneParams
    dpd_oo: dpdprops.DPDParams
    dpd_ii: dpdprops.DPDParams
    dpd_io: dpdprops.DPDParams
    dpd_rbco: dpdprops.DPDParams
    dpd_rbci: dpdprops.DPDParams
    mean_vel: float
    omega: float
    ABF_length: float
    ABF_radius: float
    B_magn: float
    m_magn: float
    CTC_diameter: float
    length_scale: pint.Quantity
    time_scale: pint.Quantity
    mass_scale: pint.Quantity

def rescale_mesh(mesh, RA: float):
    RA_current = (mesh.area / (4 * np.pi)) ** (1/2)
    scale = RA / RA_current
    mesh.vertices *= scale

def create_parameters(ureg,
                      L_,
                      R_,
                      mesh_ini,
                      mesh_ref,
                      mesh_ABF,
                      mean_vel_,
                      freq_,
                      ABF_radius_,
                      no_visc: bool=True,
                      RA: float = 5,
                      mu: float = 20,
                      Re: float = 0.5,
                      verbose: bool=True):

    @ureg.wraps(None, ureg.dimensionless)
    def to_sim(a):
        return a

    assert mean_vel_.check('[velocity]')
    assert L_.check('[length]')
    assert R_.check('[length]')
    assert freq_.check('[frequency]')

    rescale_mesh(mesh_ini, RA)
    rescale_mesh(mesh_ref, RA)

    rbc = dpdprops.JuelicherLimRBCDefaultParams(ureg)

    # physical constant
    CTC_diameter_ = 9 * ureg.um
    RA_ = 3.34 * ureg.um
    nuo_ = 1.004e-6 * ureg.m**2 / ureg.s
    etao_ = 1.002e-3 * ureg.Pa * ureg.s
    rho_ = 1e3 * ureg.kg / ureg.m**3
    kB_ = 1.38064852e-23 * ureg.m**2 *ureg.kg / (ureg.s**2 * ureg.K)
    T_ = ureg.Quantity(37, ureg.degC).to('kelvin')
    kBT_ = kB_ * T_

    # physical parameters
    omega_ = 2 * np.pi * freq_
    mean_shear_ = omega_

    # dimless numbers in "real" world
    Re_ = to_sim(mean_shear_ * RA_**2 / nuo_)
    Ca_ = to_sim(mean_shear_ * RA_ * etao_ / rbc.mu)
    E_ = to_sim(rbc.eta_m / (RA_ * etao_))
    kb_kBT_ = to_sim(rbc.kappab / kBT_)

    # go to simulation world
    rc = 1
    nd = 10
    m = 1
    rho = m * nd
    a_rc_kBTinv = 100

    Ca = Ca_
    E = E_
    kb_kBT = kb_kBT_
    visc_ratio = 5

    if Re is None:
        Re = Re_

    # solve for eta, mean_shear and etam to satisfy Re, Ca and E
    etao = np.sqrt(Ca * mu * rho * RA / Re)
    mean_shear = Ca * mu / (RA * etao)
    eta_m = E * RA * etao


    f = np.sqrt(Re / Re_)

    length_scale = RA_ / RA
    time_scale = mean_shear / (f * mean_shear_)
    mass_scale = (rho_ * length_scale**3 / rho).to(ureg.kg)

    force_scale = length_scale * mass_scale / time_scale**2

    R = to_sim(R_ / length_scale)
    L = to_sim(L_ / length_scale)


    assert length_scale.check('[length]')
    assert force_scale.check('[force]')
    assert mass_scale.check('[mass]')
    assert time_scale.check('[time]')

    assert abs(Re - RA**2 * mean_shear * rho / etao) < 1e-2

    # RBC parameters

    rbc_prms = rbc.get_params(mesh=mesh_ini,
                              length_scale=length_scale,
                              time_scale=time_scale,
                              mass_scale=mass_scale)

    if no_visc:
        rbc_prms.set_viscosity(0)
    else:
        rbc_prms.set_viscosity(eta_m)

    rbc_prms.ka /= 1e3 # Those are computational tricks to reduce the timestep limit
    rbc_prms.kv /= 1e3


    # fluid parameters

    kBT = rbc_prms.bending_modulus() / kb_kBT

    Cs = dpdprops.sound_speed(a=a_rc_kBTinv*kBT/rc,
                              nd=nd,
                              mass=m,
                              rc=rc,
                              kBT=kBT)

    U = RA * mean_shear
    Ma = U / Cs

    dpd_oo_prms = dpdprops.create_dpd_params_from_Re_Ma(Re=Re, Ma=Ma, U=U, L=RA,
                                                        nd=nd, rc=rc, mass=m,
                                                        a_rc_kBTinv=a_rc_kBTinv)

    nuo = dpd_oo_prms.kinematic_viscosity()
    nui = nuo * visc_ratio

    dpd_ii_prms = dpdprops.create_dpd_params_from_props(kinematic_viscosity=nui,
                                                        sound_speed=dpd_oo_prms.sound_speed(),
                                                        nd=nd, rc=rc, mass=m, a_rc_kBTinv=a_rc_kBTinv)

    dpd_io_prms = deepcopy(dpd_ii_prms)
    dpd_io_prms.gamma = 0 # keep only the repulsion force.


    # FSI params

    nd_2D_membrane = len(mesh_ini.vertices) / mesh_ini.area
    dpd_rbco_prms = dpdprops.create_fsi_dpd_params(fluid_params=dpd_oo_prms, nd_membrane=nd_2D_membrane)
    dpd_rbci_prms = dpdprops.create_fsi_dpd_params(fluid_params=dpd_ii_prms, nd_membrane=nd_2D_membrane)

    # ABF params

    ABF_radius = to_sim(ABF_radius_ / length_scale)
    freq = to_sim(freq_ * time_scale)
    omega = 2 * np.pi * freq
    B_magn = 10000 # TODO
    m_magn = 100 # TODO

    mesh_ABF.vertices *= ABF_radius / np.min(np.ptp(mesh_ABF.vertices, axis=0))

    ABF_extents = np.ptp(mesh_ABF.vertices, axis=0)
    ABF_extents_ = ABF_extents * length_scale
    ABF_length = max(ABF_extents)

    # CTC params

    CTC_diameter = to_sim(CTC_diameter_ / length_scale)

    mesh_CTC = trimesh.creation.icosphere(radius=CTC_diameter/2,
                                          subdivisions=3)

    mean_vel = to_sim(mean_vel_ / length_scale * time_scale)

    if verbose:
        print(f"length_scale = {length_scale.to(ureg.um)}")
        print(f"force_scale = {force_scale.to(ureg.pN)}")
        print(f"time_scale = {time_scale.to(ureg.s)}")
        print(f"mass_scale = {mass_scale.to(ureg.kg)}")

        print(f"mu = {mu}")
        print(f"etao = {etao}")
        print(f"eta_m = {eta_m}")

        print(f"Re = {Re_} (phys), {Re} (sim)")
        print(f"Ca = {Ca_} (phys), {Ca} (sim)")
        print(f"E = {E_} (phys), {E} (sim)")
        print(f"Ma = {Ma} (sim)")
        print(rbc_prms)

        print(f"dpd_oo = {dpd_oo_prms}")
        print(f"dpd_ii = {dpd_ii_prms}")
        print(f"dpd_io = {dpd_io_prms}")

        print(f"ABF dimensions = {ABF_extents_} ({ABF_extents})")
        print(f"omega = {omega_} ({omega})")
        print(f"mean_vel = {mean_vel_} ({mean_vel})")

        print(f"CTC diameter = {CTC_diameter_} ({CTC_diameter})")

        sys.stdout.flush()

    p = CapilarryFlowParams(rc=rc, m=m, nd=nd, RA=RA,
                            kBT=kBT, L=L, R=R,
                            mesh_ini=mesh_ini,
                            mesh_ref=mesh_ref,
                            mesh_ABF=mesh_ABF,
                            mesh_CTC=mesh_CTC,
                            rbc_params=rbc_prms,
                            dpd_oo=dpd_oo_prms,
                            dpd_ii=dpd_ii_prms,
                            dpd_io=dpd_io_prms,
                            dpd_rbco=dpd_rbco_prms,
                            dpd_rbci=dpd_rbci_prms,
                            ABF_length=ABF_length,
                            ABF_radius=ABF_radius,
                            mean_vel=mean_vel,
                            omega=omega,
                            B_magn=B_magn,
                            m_magn=m_magn,
                            CTC_diameter=CTC_diameter,
                            length_scale=length_scale,
                            time_scale=time_scale,
                            mass_scale=mass_scale)
    return p


def dump_parameters(p: 'CapilarryFlowParams',
                    filename: str):
    with open(filename, 'wb') as f:
        pickle.dump(p, f, pickle.HIGHEST_PROTOCOL)

def load_parameters(filename: str):
    with open(filename, 'rb') as f:
        return pickle.load(f)


def main(argv):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('mesh_ini', type=str, help="Initial RBC mesh.")
    parser.add_argument('mesh_ref', type=str, help="Stress-free RBC mesh.")
    parser.add_argument('mesh_ABF', type=str, help="ABF mesh.")
    parser.add_argument('--Re', type=float, default=0.5, help="Simulation Reynolds number.")
    parser.add_argument('--L', type=float, default=50, help="Length of the pipe, in um.")
    parser.add_argument('--R', type=float, default=10, help="Radius of the pipe, in um.")
    parser.add_argument('--mean-vel', type=float, default=0, help="mean blood velocity, in mm/s.")
    parser.add_argument('--freq', type=float, default=12, help="Frequency of rotation of the ABF, in Hz.")
    parser.add_argument('--ABF-radius', type=float, default=15, help="Length of the ABF, in um.")
    parser.add_argument('--out-prms', type=str, default="parameters.pkl", help="Save parameters to this file.")
    parser.add_argument('--out-ABF-mesh', type=str, default="ABF_mesh.ply", help="Save rescaled ABF mesh to this file.")
    parser.add_argument('--out-CTC-mesh', type=str, default="CTC_mesh.ply", help="Save rescaled CTC mesh to this file.")
    args = parser.parse_args(argv)

    mesh_ini = trimesh.load_mesh(args.mesh_ini, process=False)
    mesh_ref = trimesh.load_mesh(args.mesh_ref, process=False)

    mesh_ABF = trimesh.load_mesh(args.mesh_ABF)

    ureg = pint.UnitRegistry()

    mean_vel_ = args.mean_vel * ureg.mm / ureg.s

    L_ = args.L * ureg.um
    R_ = args.R * ureg.um
    p = create_parameters(ureg,
                          L_=L_,
                          R_=R_,
                          mesh_ini=mesh_ini,
                          mesh_ref=mesh_ref,
                          mesh_ABF=mesh_ABF,
                          mean_vel_ = mean_vel_,
                          freq_ = args.freq * ureg.Hz,
                          ABF_radius_ = args.ABF_radius * ureg.um,
                          Re=args.Re)

    dump_parameters(p, args.out_prms)

    p.mesh_ABF.export(args.out_ABF_mesh)
    p.mesh_CTC.export(args.out_CTC_mesh)



if __name__ == '__main__':
    main(sys.argv[1:])
