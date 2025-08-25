#!/usr/bin/env python

from copy import deepcopy
from dataclasses import dataclass
import dpdprops
import numpy as np
import pint
import trimesh
import sys
import pyvista as pv 
import torch
sys.path.insert(0, '..')

from utils import (sdf_read,
                   sdf_write,
                   force_read,
                   force_write)


@dataclass
class ContactParams:
    sigma: float
    eps: float
    max_contact_force: float

@dataclass
class Parameters:
    rc: float
    m: float
    nd: float
    RA: float
    kBT: float
    domain: tuple
    U: float
    mesh_ini: trimesh.Trimesh
    mesh_ref: trimesh.Trimesh
    mesh_abf: trimesh.Trimesh
    rbc_params: dpdprops.MembraneParams
    dpd_oo: dpdprops.DPDParams
    dpd_ii: dpdprops.DPDParams
    dpd_io: dpdprops.DPDParams
    dpd_rbco: dpdprops.DPDParams
    dpd_rbci: dpdprops.DPDParams
    omega: float
    abf_length: float
    abf_radius: float
    magn_B: float
    magn_m: tuple
    contact_rbcs: ContactParams
    length_scale_: pint.Quantity
    time_scale_: pint.Quantity
    mass_scale_: pint.Quantity
    magn_scale_: pint.Quantity
    shear_scale_factor: float
    wall_repulsion_length: float
    sdf_filename: str
    forces_filename: str
    path_filename : str


def rescale_mesh(mesh, RA: float):
    RA_current = (mesh.area / (4 * np.pi)) ** (1/2)
    scale = RA / RA_current
    mesh.vertices *= scale


def estimated_relative_viscosity(D: '[length]',
                                 Ht: float):
    """
    Estimate the relative viscosity from the fit in Pries 92.

    Arguments:
        D: the diameter of the pipe
        Ht: hematocrit, in [0,1]
    """
    D = D.to('um').magnitude
    eta_rel_45 = 220 * np.exp(-1.3 * D) + 3.2 - 2.44 * np.exp(-0.06 * D**0.645)
    C = (0.8 + np.exp(-0.075 * D)) * (-1 + 1 / (1 + D * 1e-11))
    eta_rel = 1 + (eta_rel_45 - 1) * ((1 - Ht)**C - 1) / ((1-0.45)**C - 1)
    return eta_rel


def create_parameters(*,
                      ureg,
                      mesh_ini: trimesh.Trimesh,
                      mesh_ref: trimesh.Trimesh,
                      mesh_abf: trimesh.Trimesh,
                      freq_: pint.Quantity,
                      abf_radius_: pint.Quantity,
                      Ht: float,
                      RA: float,
                      Re: float,
                      visc_ratio: float,
                      mu: float=20,
                      entry_radius_: pint.Quantity,
                      entry_velocity_: pint.Quantity,
                      wall_repulsion_length_: pint.Quantity,
                      magn_B_: pint.Quantity,
                      magn_m_: pint.Quantity,
                      sdf_filename: str,
                      forces_filename: str,
                      path_filename : str,
                      verbose: bool=True):

    assert freq_.check('[frequency]')

    rescale_mesh(mesh_ini, RA)
    rescale_mesh(mesh_ref, RA)

    rbc = dpdprops.JuelicherLimRBCDefaultParams(ureg)

    # physical constants
    RA_ = 3.34 * ureg.um
    nuo_ = 1.004e-6 * ureg.m**2 / ureg.s
    etao_ = 1.002e-3 * ureg.Pa * ureg.s
    rho_ = 1e3 * ureg.kg / ureg.m**3
    kB_ = 1.38064852e-23 * ureg.m**2 *ureg.kg / (ureg.s**2 * ureg.K)
    T_ = ureg.Quantity(37, ureg.degC).to('kelvin')
    kBT_ = kB_ * T_

    # physical parameters
    omega_ = 2 * np.pi * freq_
    mean_shear_ = (entry_velocity_ / entry_radius_).to('Hz')

    # dimless numbers in "real" world
    Re_ = float(mean_shear_ * RA_**2 / nuo_)
    Ca_ = float(mean_shear_ * RA_ * etao_ / rbc.mu)
    E_ = float(rbc.eta_m / (RA_ * etao_))
    kb_kBT_ = float(rbc.kappab / kBT_)
    torque_ratio_ = float(magn_B_ * magn_m_ / (abf_radius_**3 * etao_ * omega_))

    # go to simulation world
    rc = 1
    nd = 10
    m = 1
    rho = m * nd
    a_rc_kBTinv = 100

    Ca = Ca_
    E = E_
    kb_kBT = kb_kBT_

    # solve for eta, mean_shear and etam to satisfy Re, Ca and E
    etao = np.sqrt(Ca * mu * rho * RA / Re)
    mean_shear = Ca * mu / (RA * etao)
    eta_m = E * RA * etao


    f = np.sqrt(Re / Re_)

    length_scale_ = RA_ / RA
    time_scale_ = (mean_shear / (f * mean_shear_)).to('s')
    mass_scale_ = (rho_ * length_scale_**3 / rho).to(ureg.kg)

    force_scale_ = length_scale_ * mass_scale_ / time_scale_**2

    entry_radius = float(entry_radius_ / length_scale_)
    Vmax = entry_radius * mean_shear


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



    # abf params

    abf_radius = float(abf_radius_ / length_scale_)
    freq = float(freq_ * time_scale_)
    omega = 2 * np.pi * freq

    mesh_abf.vertices *= abf_radius / np.min(0.5 * np.ptp(mesh_abf.vertices, axis=0))

    abf_extents = np.ptp(mesh_abf.vertices, axis=0)
    abf_extents_ = abf_extents * length_scale_
    abf_length = max(abf_extents)

    # align the mesh to its principal axes
    I = mesh_abf.moment_inertia

    Lambda, W = trimesh.inertia.principal_axis(I)
    vertices = np.array(mesh_abf.vertices)
    vertices -= mesh_abf.center_mass
    vertices = np.dot(vertices, W.T)
    mesh_abf.vertices = vertices
    if mesh_abf.volume < 0:
        mesh_abf.faces = mesh_abf.faces[:,::-1]

    # magnetic params

    magn_B = 10        # sets the magnetic scale
    torque_ratio = torque_ratio_
    magn_m = torque_ratio * omega * abf_radius**3 * dpd_oo_prms.dynamic_viscosity() / magn_B
    magn_m = np.array([0, magn_m, 0])
    magn_m = np.dot(magn_m, W.T)
    if magn_m[2] < 0:
        magn_m *= -1
    magn_m = tuple(magn_m)

    magn_scale_ = (magn_B_ / magn_B).to(ureg.T)


    # body force
    D_ = 2 * entry_radius_
    D = 2 * entry_radius
    eta0 = dpd_oo_prms.dynamic_viscosity()
    eta_eff = eta0 * estimated_relative_viscosity(D_, Ht)
    Vmean = Vmax/2
    Q = Vmean * np.pi * entry_radius**2
    pressure_gradient = 128 * eta_eff * Q / (np.pi * D**4)


    # SDF
    scale_path= 0.041784223
    path_numpy = np.load(path_filename)/scale_path
    
    n, h, sdf = sdf_read(sdf_filename)
    max_rad = abs(np.min(sdf))
    scale = entry_radius / max_rad
    sdf *= scale
    path = scale *path_numpy
    rescaled_path_filename = 'path_rescaled.npy'
    np.save(rescaled_path_filename,path)
    points = path.astype(float)
    line = pv.Spline(points, n_points=len(points))
    line.save('paraview_path.vtp')

    h = [scale * h_ for h_ in h]
    rescaled_sdf_filename="retina.rescaled.sdf"
    sdf_write(rescaled_sdf_filename, sdf, n, h)

    domain=(n[0] * h[0],
            n[1] * h[1],
            n[2] * h[2])

    # forces
    n, h, forces = force_read(forces_filename)
    magn = np.linalg.norm(forces, axis=-1)
    fscale = pressure_gradient / np.max(magn) / nd
    fscale *= -4 # TODO: factor to be tuned
    forces *= fscale
    h = [scale * h_ for h_ in h]
    rescaled_force_filename="forces.rescaled.sdf"
    force_write(rescaled_force_filename, forces, n, h)


    # Contact forces

    max_contact_force = dpd_oo_prms.a * 50

    l0 = rbc_prms.get_l0()
    sigma = 4 * l0 / 3
    eps = 50 * rbc_prms.bending_modulus()
    contact_rbcs = ContactParams(sigma=sigma,
                                 eps=eps,
                                 max_contact_force=max_contact_force)

    wall_repulsion_length = float(wall_repulsion_length_ / length_scale_)

    if verbose:
        print(f"length_scale = {length_scale_.to(ureg.um)}")
        print(f"force_scale = {force_scale_.to(ureg.pN)}")
        print(f"time_scale = {time_scale_.to(ureg.s)}")
        print(f"mass_scale = {mass_scale_.to(ureg.kg)}")
        print(f"magn_scale = {magn_scale_.to(ureg.T)}")
        print()

        print(f"Domain size: {domain}")
        print()

        print(f"mu = {mu}")
        print(f"Ventry = {entry_velocity_} ({float(entry_velocity_/length_scale_*time_scale_)})")
        print(f"etao = {etao_} ({etao})")
        print(f"mean shear = {mean_shear_.to(1/ureg.s)} ({mean_shear})")
        print(f"force field scale = {fscale}")
        print()

        print(f"Re = {Re_} (phys), {Re} (sim)")
        print(f"Ca = {Ca_} (phys), {Ca} (sim)")
        print(f"E = {E_} (phys), {E} (sim)")
        print(f"Ma = {Ma} (sim)")
        print(rbc_prms)

        print(f"dpd_oo = {dpd_oo_prms}")
        print(f"dpd_ii = {dpd_ii_prms}")
        print(f"dpd_io = {dpd_io_prms}")
        print(f"dpd_rbco = {dpd_rbco_prms}")
        print(f"dpd_rbci = {dpd_rbci_prms}")

        print(f"pressure_gradient = {pressure_gradient}")
        print(f"contact params = {contact_rbcs}")
        print()

        print(f"abf dimensions = {abf_extents_} ({abf_extents})")
        print(f"omega = {omega_} ({omega})")
        print(f"B = {magn_B_} ({magn_B})")
        print(f"m = {magn_m_} ({magn_m})")

        sys.stdout.flush()

    p = Parameters(rc=rc, m=m, nd=nd, RA=RA,
                   kBT=kBT, U=U,
                   domain=domain,
                   mesh_ini=mesh_ini,
                   mesh_ref=mesh_ref,
                   mesh_abf=mesh_abf,
                   rbc_params=rbc_prms,
                   dpd_oo=dpd_oo_prms,
                   dpd_ii=dpd_ii_prms,
                   dpd_io=dpd_io_prms,
                   dpd_rbco=dpd_rbco_prms,
                   dpd_rbci=dpd_rbci_prms,
                   abf_length=abf_length,
                   abf_radius=abf_radius,
                   omega=omega,
                   magn_B=magn_B,
                   magn_m=magn_m,
                   contact_rbcs=contact_rbcs,
                   length_scale_=length_scale_,
                   time_scale_=time_scale_,
                   mass_scale_=mass_scale_,
                   magn_scale_=magn_scale_,
                   shear_scale_factor=f,
                   wall_repulsion_length=wall_repulsion_length,
                   sdf_filename=rescaled_sdf_filename,
                   forces_filename=rescaled_force_filename,
                   path_filename=rescaled_path_filename)
    return p


def dump_parameters(p: 'Parameters',
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
    parser.add_argument('mesh_abf', type=str, help="abf mesh.")
    parser.add_argument('--rbc-res', type=int, default=3, help="Subdivision level of the mesh.")
    parser.add_argument('--Re', type=float, default=0.5, help="Simulation Reynolds number.")
    parser.add_argument('--RA', type=float, default=4, help="Reduced radius of the RBC, in simulation units.")
    parser.add_argument('--Ht', type=float, default=0.05, help="Hematocrit")
    parser.add_argument('--visc-ratio', type=float, default=5, help="Viscosity ratio between inner and outer solvent")
    parser.add_argument('--wall-repulsion-length', type=str, default="0.5_um", help="Repulsion length of the wall forces")
    parser.add_argument('--sdf', type=str, required=True, help="SDF file")
    parser.add_argument('--forces', type=str, required=True, help="forces file")
    parser.add_argument('--entry-radius', type=str, default='5_um', help="radius of largest pipe")
    parser.add_argument('--entry-velocity', type=str, default='1_mm_per_s', help="inlet velocity")
    parser.add_argument('--freq', type=str, default="10_Hz", help="Frequency of rotation of the abf.")
    parser.add_argument('--abf-radius', type=str, default="2_um", help="Radius of the abf.")
    parser.add_argument('--magn-m', type=str, default="1e-14_N_m_per_T", help="Magnitude of the magnetic moment of the swimmer.")
    parser.add_argument('--out-params', type=str, default="parameters.pkl", help="Save all simulation parameters to this file.")
    parser.add_argument('--out-abf-mesh', type=str, default="abf_mesh.ply", help="Save rescaled abf mesh to this file.")
    parser.add_argument('--dir_path', type=str, required=True, help="Path A* file")
    args = parser.parse_args(argv)

    rbc_res = args.rbc_res
    mesh_ini = dpdprops.load_equilibrium_mesh(rbc_res)
    mesh_ref = dpdprops.load_stress_free_mesh(rbc_res)

    mesh_abf = trimesh.load_mesh(args.mesh_abf)

    ureg = pint.UnitRegistry()

    wall_repulsion_length_ = get_quantity(ureg, args.wall_repulsion_length)
    entry_radius_ = get_quantity(ureg, args.entry_radius)
    entry_velocity_ = get_quantity(ureg, args.entry_velocity)

    abf_radius_ = get_quantity(ureg, args.abf_radius)

    freq_ = get_quantity(ureg, args.freq)
    magn_B_ = 3 * ureg.mT
    magn_m_ = get_quantity(ureg, args.magn_m)

    p = create_parameters(ureg=ureg,
                          Re=args.Re,
                          RA=args.RA,
                          Ht=args.Ht,
                          visc_ratio=args.visc_ratio,
                          mesh_ini=mesh_ini,
                          mesh_ref=mesh_ref,
                          mesh_abf=mesh_abf,
                          freq_=freq_,
                          magn_B_=magn_B_,
                          magn_m_=magn_m_,
                          abf_radius_=abf_radius_,
                          wall_repulsion_length_=wall_repulsion_length_,
                          entry_velocity_=entry_velocity_,
                          entry_radius_=entry_radius_,
                          sdf_filename=args.sdf,
                          forces_filename=args.forces,
                          path_filename = args.dir_path,
                          )
    
    dump_parameters(p, args.out_params)
    p.mesh_abf.export(args.out_abf_mesh)


if __name__ == '__main__':
    main(sys.argv[1:])
