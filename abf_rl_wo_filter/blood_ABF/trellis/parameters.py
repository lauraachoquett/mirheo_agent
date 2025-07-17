#!/usr/bin/env python

from copy import deepcopy
from dataclasses import dataclass
import numpy as np
import pickle
import pint
import trimesh
import sys

import dpdprops

from geometry import create_grid

@dataclass
class SimulationParams:
    rc: float
    m: float
    nd: float
    RA: float
    kBT: float
    domain: tuple
    L: float
    R: float
    num_bif_x: int
    num_bif_y: int
    theta_deg: float
    mesh_ini: trimesh.Trimesh
    mesh_ref: trimesh.Trimesh
    mesh_ABF: trimesh.Trimesh
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
    magn_B: float
    magn_m: tuple
    length_scale_: pint.Quantity
    time_scale_: pint.Quantity
    mass_scale_: pint.Quantity


def rescale_mesh(mesh, RA: float):
    RA_current = (mesh.area / (4 * np.pi)) ** (1/2)
    scale = RA / RA_current
    mesh.vertices *= scale


def create_parameters(ureg,
                      L_: pint.Quantity,
                      R_: pint.Quantity,
                      num_bif_x: int,
                      num_bif_y: int,
                      theta_deg: float,
                      mesh_ini: trimesh.Trimesh,
                      mesh_ref: trimesh.Trimesh,
                      mesh_ABF: trimesh.Trimesh,
                      mean_vel_: pint.Quantity,
                      freq_: pint.Quantity,
                      magn_B_: pint.Quantity,
                      magn_m_: pint.Quantity,
                      ABF_radius_: pint.Quantity,
                      RA: float,
                      Re: float = 0.5,
                      mu: float = 20,
                      verbose: bool=True):

    assert mean_vel_.check('[velocity]')
    assert L_.check('[length]')
    assert R_.check('[length]')
    assert freq_.check('[frequency]')

    rescale_mesh(mesh_ini, RA)
    rescale_mesh(mesh_ref, RA)

    rbc = dpdprops.JuelicherLimRBCDefaultParams(ureg)

    # physical constants
    RA_ = dpdprops.equivalent_sphere_radius(area=rbc.A0)
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
    Re_ = float(mean_shear_ * RA_**2 / nuo_)
    Ca_ = float(mean_shear_ * RA_ * etao_ / rbc.mu)
    E_ = float(rbc.eta_m / (RA_ * etao_))
    kb_kBT_ = float(rbc.kappab / kBT_)
    torque_ratio_ = float(magn_B_ * magn_m_ / (ABF_radius_**3 * etao_ * omega_))

    # go to simulation world
    rc = 1
    nd = 10
    m = 1
    rho = m * nd
    a_rc_kBTinv = 100

    Ca = Ca_
    E = E_
    kb_kBT = kb_kBT_
    torque_ratio = torque_ratio_
    visc_ratio = 5

    if Re is None:
        Re = Re_

    # solve for eta, mean_shear and etam to satisfy Re, Ca and E
    etao = np.sqrt(Ca * mu * rho * RA / Re)
    mean_shear = Ca * mu / (RA * etao)
    eta_m = E * RA * etao


    f = np.sqrt(Re / Re_)

    length_scale_ = RA_ / RA
    time_scale_ = mean_shear / (f * mean_shear_)
    mass_scale_ = (rho_ * length_scale_**3 / rho).to(ureg.kg)

    force_scale_ = length_scale_ * mass_scale_ / time_scale_**2

    L = float(L_ / length_scale_)
    R = float(R_ / length_scale_)


    assert length_scale_.check('[length]')
    assert force_scale_.check('[force]')
    assert mass_scale_.check('[mass]')
    assert time_scale_.check('[time]')

    assert abs(Re - RA**2 * mean_shear * rho / etao) < 1e-2

    # RBC parameters

    rbc_prms = rbc.get_params(mesh=mesh_ini,
                              length_scale=length_scale_,
                              time_scale=time_scale_,
                              mass_scale=mass_scale_)

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

    ABF_radius = float(ABF_radius_ / length_scale_)
    freq = float(freq_ * time_scale_)
    omega = 2 * np.pi * freq

    mesh_ABF.vertices *= 2 * ABF_radius / np.min(np.ptp(mesh_ABF.vertices, axis=0))

    ABF_extents = np.ptp(mesh_ABF.vertices, axis=0)
    ABF_extents_ = ABF_extents * length_scale_
    ABF_length = max(ABF_extents)

    mean_vel = float(mean_vel_ / length_scale_ * time_scale_)

    # align the mesh to its principal axes
    I = mesh_ABF.moment_inertia

    Lambda, W = trimesh.inertia.principal_axis(I)
    vertices = np.array(mesh_ABF.vertices)
    vertices -= mesh_ABF.center_mass
    vertices = np.dot(vertices, W.T)
    mesh_ABF.vertices = vertices
    if mesh_ABF.volume < 0:
        mesh_ABF.faces = mesh_ABF.faces[:,::-1]


    # magnetic params

    magn_B = 10        # sets the magnetic scale
    magn_m = torque_ratio * omega * ABF_radius**3 * dpd_oo_prms.dynamic_viscosity() / magn_B
    magn_m = np.array([0, magn_m, 0])
    magn_m = np.dot(magn_m, W.T)
    if magn_m[2] < 0:
        magn_m *= -1
    magn_m = tuple(magn_m)

    magn_scale_ = (magn_B_ / magn_B).to(ureg.T)


    # Geometry

    sdf_grid, domain = create_grid(L=L, R=R,
                                   nx=num_bif_x,
                                   ny=num_bif_y,
                                   theta_deg=theta_deg,
                                   margin=2*rc,
                                   h=0.25*rc)

    if verbose:
        print(f"domain: {domain}")
        print(f"length_scale = {length_scale_.to(ureg.um)}")
        print(f"force_scale = {force_scale_.to(ureg.pN)}")
        print(f"time_scale = {time_scale_.to(ureg.s)}")
        print(f"mass_scale = {mass_scale_.to(ureg.kg)}")

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

        sys.stdout.flush()

    p = SimulationParams(rc=rc, m=m, nd=nd, RA=RA,
                         domain=domain,
                         kBT=kBT,
                         L=L, R=R,
                         num_bif_x=num_bif_x,
                         num_bif_y=num_bif_y,
                         theta_deg=theta_deg,
                         mesh_ini=mesh_ini,
                         mesh_ref=mesh_ref,
                         mesh_ABF=mesh_ABF,
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
                         magn_B=magn_B,
                         magn_m=magn_m,
                         length_scale_=length_scale_,
                         time_scale_=time_scale_,
                         mass_scale_=mass_scale_)
    return sdf_grid, p


def dump_parameters(p: 'SimulationParams',
                    filename: str):
    with open(filename, 'wb') as f:
        pickle.dump(p, f, pickle.HIGHEST_PROTOCOL)


def load_parameters(filename: str):
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
    parser.add_argument('mesh_ABF', type=str, help="ABF mesh.")
    parser.add_argument('--rbc_mesh_res', type=int, default=3, help="Subdivision level of the RBC mesh.")
    parser.add_argument('--Re', type=float, default=0.5, help="Simulation Reynolds number.")
    parser.add_argument('--RA', type=float, default=5, help="Reduced radius of the RBC in simulation units.")
    parser.add_argument('--L',  type=str, default="30_um", help="Pipe length between two bifurcations.")
    parser.add_argument('--R',  type=str, default="6_um", help="Radius of the main pipe.")
    parser.add_argument('--num-bif-x',  type=int, default=3, help="Number of bifurcations along x.")
    parser.add_argument('--num-bif-y',  type=int, default=3, help="Number of bifurcations along y.")
    parser.add_argument('--theta-deg',  type=float, default=60, help="Angle between the branches of the bifurcations, in degrees.")
    parser.add_argument('--mean-vel',   type=str, default="0_mm_per_s", help="mean blood velocity.")
    parser.add_argument('--freq',       type=str, default="100_Hz", help="Frequency of rotation of the ABF.")
    parser.add_argument('--magn-m', type=str, default="1e-14_N_m_per_T", help="Magnitude of the magnetic moment of the swimmer.")
    parser.add_argument('--ABF-radius', type=str, default="2_um", help="Radius of the ABF.")
    parser.add_argument('--out-prms',   type=str, default="parameters.pkl", help="Save parameters to this file.")
    parser.add_argument('--out-ABF-mesh', type=str, default="ABF_mesh.ply", help="Save rescaled ABF mesh to this file.")
    parser.add_argument('--out-geometry', type=str, nargs='+', default=["trellis.sdf"], help="Save the SDF of the geometry to these files.")
    args = parser.parse_args(argv)

    rbc_res = args.rbc_mesh_res
    mesh_ini = dpdprops.load_equilibrium_mesh(rbc_res)
    mesh_ref = dpdprops.load_stress_free_mesh(rbc_res)

    mesh_ABF = trimesh.load_mesh(args.mesh_ABF)

    ureg = pint.UnitRegistry()

    mean_vel_ = get_quantity(ureg, args.mean_vel)
    L_ = get_quantity(ureg, args.L)
    R_ = get_quantity(ureg, args.R)
    freq_ = get_quantity(ureg, args.freq)
    ABF_radius_ = get_quantity(ureg, args.ABF_radius)

    magn_B_ = 3 * ureg.mT
    magn_m_ = get_quantity(ureg, args.magn_m)

    sdf_grid, p = create_parameters(ureg,
                                    L_=L_, R_=R_,
                                    num_bif_x=args.num_bif_x,
                                    num_bif_y=args.num_bif_y,
                                    theta_deg=args.theta_deg,
                                    mesh_ini=mesh_ini,
                                    mesh_ref=mesh_ref,
                                    mesh_ABF=mesh_ABF,
                                    mean_vel_=mean_vel_,
                                    freq_=freq_,
                                    magn_B_=magn_B_,
                                    magn_m_=magn_m_,
                                    ABF_radius_=ABF_radius_,
                                    Re=args.Re,
                                    RA=args.RA)

    dump_parameters(p, args.out_prms)

    p.mesh_ABF.export(args.out_ABF_mesh)

    for f in args.out_geometry:
        sdf_grid.dump(f)



if __name__ == '__main__':
    main(sys.argv[1:])
