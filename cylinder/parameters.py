from copy import deepcopy
from dataclasses import dataclass
import mirheo as mir
from dpdprops import (
    JuelicherLimRBCDefaultParams,
    create_dpd_params_from_props,
    create_fsi_dpd_params,
    load_stress_free_mesh,
    load_equilibrium_mesh,
    sound_speed,
    create_dpd_params_from_Re_Ma,
)
import pint
import trimesh
import numpy as np

import argparse


@dataclass
class ContactParams:
    sigma: float
    eps: float
    max_contact_force: float


def rescale_by_area(mesh, A):
    mesh.vertices *= np.sqrt(A / mesh.area)
    return mesh


def get_quantity(ureg: pint.UnitRegistry, s: str):
    if "/" in s or " " in s:
        raise ValueError(
            f"'{s}' should not contain '/' or spaces, otherwise it will create directories with undesired names. Use '_per_' and '_' instead."
        )
    return ureg(s.replace("_", " "))


def estimated_relative_viscosity(D: "[length]", Ht: float):
    """
    Estimate the relative viscosity from the fit in Pries 92.

    Arguments:
        D: the diameter of the pipe
        Ht: hematocrit, in [0,1]
    """
    D = D.to("um").magnitude
    eta_rel_45 = 220 * np.exp(-1.3 * D) + 3.2 - 2.44 * np.exp(-0.06 * D**0.645)
    C = (0.8 + np.exp(-0.075 * D)) * (-1 + 1 / (1 + D * 1e-11))
    eta_rel = 1 + (eta_rel_45 - 1) * ((1 - Ht) ** C - 1) / ((1 - 0.45) ** C - 1)
    return eta_rel


def set_parameters(Re, verbose):
    parser = argparse.ArgumentParser(
        description="Run RBCs flowing in a capillary pipe."
    )
    parser.add_argument("--rbc-res", type=int, help="Subdivision level of the mesh.")
    parser.add_argument("--L", type=str, default="50_um", help="Length of the pipe.")
    parser.add_argument("--R", type=str, default="10_um", help="Radius of the pipe.")
    parser.add_argument(
        "--Re", type=float, default=Re, help="Simulation Reynolds number."
    )
    parser.add_argument(
        "--RA",
        type=float,
        default=6,
        help="Reduced radius of the RBC, in simulation units.",
    )
    parser.add_argument(
        "--Vmax", type=str, default="6.74_mm_per_s", help="Maximum flow velocity."
    )
    parser.add_argument("--Ht", type=float, default=0.15, help="Hematocrit")
    parser.add_argument(
        "--C",
        type=float,
        default=5,
        help="Viscosity ratio between inner and outer solvent",
    )
    parser.add_argument(
        "--wall-repulsion-length",
        type=str,
        default="0.5_um",
        help="Repulsion length of the wall forces",
    )
    parser.add_argument(
        "--out-params",
        type=str,
        default="parameters.pkl",
        help="Save all simulation parameters to this file.",
    )
    args = parser.parse_args()

    ureg = pint.UnitRegistry()

    Vmax_ = get_quantity(ureg, args.Vmax)
    L_ = get_quantity(ureg, args.L)
    R_ = get_quantity(ureg, args.R)

    wall_repulsion_length_ = get_quantity(ureg, args.wall_repulsion_length)

    params_juelicher = JuelicherLimRBCDefaultParams(ureg)
    V0 = params_juelicher.V0

    # physical constants
    RA_ = 3.34 * ureg.um
    nuo_ = 1.004e-6 * ureg.m**2 / ureg.s
    etao_ = 1.002e-3 * ureg.Pa * ureg.s
    rho_ = 1e3 * ureg.kg / ureg.m**3
    kB_ = 1.38064852e-23 * ureg.m**2 * ureg.kg / (ureg.s**2 * ureg.K)
    T_ = ureg.Quantity(37, ureg.degC).to("kelvin")
    kBT_ = kB_ * T_

    # physical parameters
    mean_shear_ = Vmax_ / R_
    # dimless numbers in "real" world
    Re_ = float(mean_shear_ * RA_**2 / nuo_)
    Ca_ = float(mean_shear_ * RA_ * etao_ / params_juelicher.mu)
    E_ = float(params_juelicher.eta_m / (RA_ * etao_))
    kb_kBT_ = float(params_juelicher.kappab / kBT_)

    # go to simulation world
    rc = 1
    nd = 10
    Re = args.Re
    RA = args.RA
    Ht = args.Ht
    visc_ratio = args.C

    mass = 1
    rho = mass * nd
    mu = 20
    a_rc_kBTinv = 100
    Ca = Ca_
    E = E_
    kb_kBT = kb_kBT_

    # solve for eta, mean_shear and etam to satisfy Re, Ca and E
    etao = np.sqrt(Ca * mu * rho * RA / Re)
    mean_shear = Ca * mu / (RA * etao)
    eta_m = E * RA * etao

    # Scale :
    f = np.sqrt(Re / Re_)
    length_scale_ = RA_ / RA
    time_scale_ = (mean_shear / (f * mean_shear_)).to("s")
    mass_scale_ = (rho_ * length_scale_**3 / rho).to(ureg.kg)
    force_scale_ = length_scale_ * mass_scale_ / time_scale_**2

    # Cylinder
    R = float(R_ / length_scale_)
    L = float(L_ / length_scale_)
    Vmax = R * mean_shear

    # Mesh RBC + INIT
    mesh_ini = load_equilibrium_mesh()
    mesh_ref = load_stress_free_mesh()

    A0 = 4 * np.pi * RA**2
    mesh_ini = rescale_by_area(mesh_ini, A0)
    mesh_ref = rescale_by_area(mesh_ref, A0)

    rbc_params = params_juelicher.get_params(
        length_scale=length_scale_,
        time_scale=time_scale_,
        mass_scale=mass_scale_,
        mesh=mesh_ini,
    )

    rbc_params.set_viscosity(eta_m)
    rbc_params.ka /= 1e3  # Those are computational tricks to reduce the timestep limit
    rbc_params.kv /= 1e3

    # Fluid
    kBT = rbc_params.bending_modulus() / kb_kBT

    Cs = sound_speed(a=a_rc_kBTinv * kBT / rc, nd=nd, mass=mass, rc=rc, kBT=kBT)

    U = RA * mean_shear
    Ma = U / Cs

    params_outer = create_dpd_params_from_Re_Ma(
        Re=Re, Ma=Ma, U=U, L=RA, nd=nd, rc=rc, mass=mass, a_rc_kBTinv=a_rc_kBTinv
    )
    nuo = params_outer.kinematic_viscosity()
    nui = nuo * visc_ratio

    params_inner = create_dpd_params_from_props(
        kinematic_viscosity=nui,
        sound_speed=params_outer.sound_speed(),
        nd=nd,
        rc=rc,
        mass=mass,
        a_rc_kBTinv=a_rc_kBTinv,
    )

    params_inner_outer = deepcopy(params_inner)
    params_inner_outer.gamma = 0

    # FSI
    nd_2D_membrane = len(mesh_ini.vertices) / mesh_ini.area

    inner_rbc_params = create_fsi_dpd_params(
        fluid_params=params_inner, nd_membrane=nd_2D_membrane
    )
    outer_rbc_params = create_fsi_dpd_params(
        fluid_params=params_outer, nd_membrane=nd_2D_membrane
    )

    ## Timestep :
    dt_fluid = min(
        params_outer.get_max_dt(),
        params_inner.get_max_dt(),
        inner_rbc_params.get_max_dt(),
        outer_rbc_params.get_max_dt(),
    )

    dt_rbc_el = rbc_params.get_max_dt_elastic(mass=mass)
    dt_rbc_visc = rbc_params.get_max_dt_visc(mass=mass)
    dt_rbc = rbc_params.get_max_dt(mass=mass)
    substeps = 5 + int(dt_fluid / dt_rbc_el)
    dt = dt_fluid

    # body force
    D_ = 2 * R_
    D = 2 * R
    eta0 = params_outer.dynamic_viscosity()
    eta_eff = eta0 * estimated_relative_viscosity(D_, Ht)
    Vmean = Vmax / 2
    Q = Vmean * np.pi * R**2
    pressure_gradient = 128 * eta_eff * Q / (np.pi * D**4)
    f = Vmax * 4 * eta_eff / (R**2 * nd)

    # Contact forces
    l0 = rbc_params.get_l0()
    sigma = 4 * l0 / 3
    eps = 50 * rbc_params.bending_modulus()
    max_contact_force = params_outer.a * 50
    contact_rbcs = ContactParams(
        sigma=sigma, eps=eps, max_contact_force=max_contact_force
    )

    if verbose:
        print("\n--- Physical Data ---")
        print(f"RA_ (physical RBC radius): {RA_}")
        print(f"nuo_ (physical kinematic viscosity): {nuo_}")
        print(f"etao_ (physical dynamic viscosity): {etao_}")
        print(f"rho_ (physical density): {rho_}")
        print(f"kB_ (Boltzmann constant): {kB_}")
        print(f"T_ (temperature): {T_}")
        print(f"kBT_ (thermal energy): {kBT_}")
        print(f"mean_shear_ (physical mean shear): {mean_shear_}")
        print(f"Re_ (physical Reynolds number): {Re_}")
        print(f"Ca_ (physical Capillary number): {Ca_}")
        print(f"E_ (physical E): {E_}")
        print(f"kb_kBT_ (physical kb/kBT): {kb_kBT_}")
        print(f"L_ (pipe length): {L_}")
        print(f"R_ (pipe radius): {R_}")
        print(f"Vmax_ (max velocity): {Vmax_}")

        print("\n--- Simulation Data ---")
        print(f"rc (DPD cutoff): {rc}")
        print(f"nd (DPD number density): {nd}")
        print(f"Re (simulation Reynolds number): {Re}")
        print(f"RA (simulation RBC radius): {RA}")
        print(f"Ht (hematocrit): {Ht}")
        print(f"visc_ratio (viscosity ratio): {visc_ratio}")
        print(f"mass: {mass}")
        print(f"rho: {rho}")
        print(f"mu: {mu}")
        print(f"a_rc_kBTinv: {a_rc_kBTinv}")
        print(f"E: {E}")
        print(f"Ca: {Ca}")
        print(f"Ma: {Ma}")
        print(f"kb_kBT: {kb_kBT}")
        print(f"etao (sim dynamic viscosity): {etao}")
        print(f"mean_shear (sim mean shear): {mean_shear}")
        print(f"eta_m (membrane viscosity): {eta_m}")
        print(f"f (scaling factor): {f}")
        print(f"length_scale_: {length_scale_}")
        print(f"time_scale_: {time_scale_}")
        print(f"mass_scale_: {mass_scale_}")
        print(f"force_scale_: {force_scale_}")
        print(f"R (sim pipe radius): {R}")
        print(f"L (sim pipe length): {L}")
        print(f"Vmax (sim max velocity): {Vmax}")
        print(f"dt_rbc: {dt_rbc}")
        print(f"dt_dpd: {dt_fluid}")
        print(f"dt0: {dt}")
        print(f"eta0 (sim dynamic viscosity): {eta0}")
        print(f"eta_eff (effective viscosity): {eta_eff}")
        print(f"Vmean: {Vmean}")
        print(f"Q (flow rate): {Q}")
        print(f"pressure_gradient: {pressure_gradient}")
        print(f"body force f: {f}")
        print(f"Contact sigma: {sigma}")
        print(f"Contact eps: {eps}")
        print(f"Contact max_contact_force: {max_contact_force}")

    return {
        "rbc_params": rbc_params,
        "params_outer": params_outer,
        "params_inner": params_inner,
        "params_inner_outer": params_inner_outer,
        "inner_rbc_params": inner_rbc_params,
        "outer_rbc_params": outer_rbc_params,
        "contact_rbcs": contact_rbcs,
        "dt0": dt,
        "substeps": substeps,
        "mesh_ini": mesh_ini,
        "mesh_ref": mesh_ref,
        "nd": nd,
        "mass": mass,
        "rc": rc,
        "V0": rbc_params.volume,
        "RA": RA,
        "R": R,
        "l": L,
        "f": f,
        "Vmax": Vmax,
    }
