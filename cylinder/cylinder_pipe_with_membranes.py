#!/usr/bin/env python

import mirheo as mir
from dpdprops import create_dpd_params_from_Re_Ma
import numpy as np
import os
import h5py as h5
import glob
import matplotlib.pyplot as plt
from datetime import datetime
import json
from dpdprops import KantorWLCRBCDefaultParams
import pint
import trimesh
from math import sqrt, ceil
from parameters import set_parameters

ureg = pint.UnitRegistry()


def cylinder_flow(ranks, path, params, t_end):

    rbc_params = params["rbc_params"]
    params_outer = params["params_outer"]
    params_inner = params["params_inner"]
    params_inner_outer = params["params_inner_outer"]
    inner_rbc_params = params["inner_rbc_params"]
    outer_rbc_params = params["outer_rbc_params"]
    contact_rbcs = params["contact_rbcs"]

    dt = params["dt0"]
    substeps = params["substeps"]

    mesh_ini = params["mesh_ini"]
    mesh_ref = params["mesh_ref"]

    nd = params["nd"]
    mass = params["mass"]
    rc = params["rc"]
    V0 = params["V0"]
    RA = params["RA"]
    f = params["f"]

    R = params["R"]
    l = params["l"]

    sigma = contact_rbcs.sigma
    eps = contact_rbcs.eps
    max_contact_force = contact_rbcs.max_contact_force

    print("Simulation end time : ", t_end * dt)

    print("Compute Domain and list of RBC centers")
    list_center, domain, radius_cy, center_cy = parameters_cylinder(
        V0, rc, RA, 0.15, l, R
    )
    print("Domain :", domain)
    ## Create path
    os.makedirs(path, exist_ok=True)

    u = mir.Mirheo(ranks, domain, debug_level=0, log_filename="log")
    mesh_rbc = mir.ParticleVectors.MembraneMesh(
        mesh_ini.vertices.tolist(), mesh_ref.vertices.tolist(), mesh_ini.faces.tolist()
    )

    pv_rbc = mir.ParticleVectors.MembraneVector("rbc", mass=mass, mesh=mesh_rbc)
    # place initial membrane
    # we need a position pos and an orientation described by a quaternion q
    # here we create only one membrane at the center of the domain

    quat = [1.0, 0.0, 1.0, 0.0]

    pos_qs = []
    for center_rbc in list_center:
        pos_q = list(center_rbc) + list(quat)
        pos_qs.append(pos_q)

    ic_rbc = mir.InitialConditions.Membrane(pos_qs)
    u.registerParticleVector(pv_rbc, ic_rbc)

    ### Create solvent :
    pv_solvent = mir.ParticleVectors.ParticleVector("pv", mass=mass)
    ic_solvent = mir.InitialConditions.Uniform(nd)
    u.registerParticleVector(pv_solvent, ic_solvent)

    ### Membrane solvent :
    inner_checker = mir.BelongingCheckers.Mesh("inner_solvent_checker")
    u.registerObjectBelongingChecker(inner_checker, pv_rbc)
    pv_inner = u.applyObjectBelongingChecker(
        inner_checker, pv_solvent, correct_every=0, inside="pv_inner"
    )

    int_rbc = mir.Interactions.MembraneForces(
        "int_rbc", **rbc_params.to_interactions_zero_visc(), stress_free=True
    )

    int_dpd_oo = mir.Interactions.Pairwise(
        "dpd_oo", rc, kind="DPD", **params_outer.to_interactions()
    )
    int_dpd_ii = mir.Interactions.Pairwise(
        "dpd_ii", rc, kind="DPD", **params_inner.to_interactions()
    )
    int_dpd_io = mir.Interactions.Pairwise(
        "dpd_io", rc, kind="DPD", **params_inner_outer.to_interactions()
    )
    int_dpd_mi = mir.Interactions.Pairwise(
        "dpd_mi", rc, kind="DPD", **inner_rbc_params.to_interactions()
    )
    int_dpd_mo = mir.Interactions.Pairwise(
        "dpd_mo", rc, kind="DPD", **outer_rbc_params.to_interactions()
    )

    cnt = mir.Interactions.Pairwise(
        "contact",
        rc,
        kind="RepulsiveLJ",
        epsilon=eps,
        sigma=sigma,
        aware_mode="Object",
        max_force=max_contact_force,
    )

    u.registerInteraction(int_rbc)

    u.registerInteraction(int_dpd_oo)
    u.registerInteraction(int_dpd_ii)
    u.registerInteraction(int_dpd_io)
    u.registerInteraction(int_dpd_mi)
    u.registerInteraction(int_dpd_mo)

    u.registerInteraction(cnt)

    vv = mir.Integrators.VelocityVerlet_withConstForce("vv", force=(f, 0, 0))

    wall = mir.Walls.Cylinder(
        "cylinder", center=center_cy, radius=radius_cy, axis="x", inside=True
    )
    u.registerWall(wall)  # register the wall in the coordinator
    pv_frozen = u.makeFrozenWallParticles(
        "frozen",
        walls=[wall],
        interactions=[int_dpd_oo],
        integrator=vv,
        number_density=nd,
        nsteps=1000,
        dt=params_outer.get_max_dt(),
    )

    u.setWall(wall, pv_solvent)

    u.setInteraction(int_dpd_oo, pv_solvent, pv_solvent)
    u.setInteraction(int_dpd_ii, pv_inner, pv_inner)
    u.setInteraction(int_dpd_io, pv_inner, pv_solvent)

    # Wall :
    u.setInteraction(int_dpd_oo, pv_solvent, pv_frozen)
    u.setInteraction(int_dpd_io, pv_inner, pv_frozen)
    u.setInteraction(int_dpd_mo, pv_rbc, pv_frozen)

    # Membrane :
    u.setInteraction(int_dpd_mo, pv_solvent, pv_rbc)
    u.setInteraction(int_dpd_mi, pv_inner, pv_rbc)

    u.setInteraction(cnt, pv_rbc, pv_rbc)
    u.setInteraction(cnt, pv_rbc, pv_frozen)

    u.registerIntegrator(vv)
    u.setIntegrator(vv, pv_solvent)
    u.setIntegrator(vv, pv_inner)

    vv_rbc = mir.Integrators.SubStep("vv_rbc", substeps=substeps, fastForces=[int_rbc])
    u.registerIntegrator(vv_rbc)

    u.setIntegrator(vv_rbc, pv_rbc)

    # Bouncers
    ##########
    bouncer = mir.Bouncers.Mesh("membrane_bounce", "bounce_maxwell", kBT=0.5)
    # we register the bouncer object as any other object
    u.registerBouncer(bouncer)

    # now we can set what PVs bounce on what OV:
    u.setBouncer(bouncer, pv_rbc, pv_solvent)
    u.setBouncer(bouncer, pv_rbc, pv_inner)

    dump_every = 10000
    u.registerPlugins(mir.Plugins.createStats("stats", every=dump_every))

    u.registerPlugins(
        mir.Plugins.createDumpParticles(
            "part_dump_inner",
            pv_inner,
            dump_every,
            [],
            os.path.join(path, "inner/inner-"),
        )
    )
    u.registerPlugins(
        mir.Plugins.createDumpParticles(
            "part_dump_solvent",
            pv_solvent,
            dump_every,
            [],
            os.path.join(path, "solvent/solvent-"),
        )
    )
    u.registerPlugins(
        mir.Plugins.createDumpMesh(
            "mesh_dump", pv_rbc, dump_every, os.path.join(path, "membranes_ply/")
        )
    )

    # we can also dump the frozen particles for visualization purpose
    u.registerPlugins(
        mir.Plugins.createDumpParticles(
            "part_dump_wall",
            pv_frozen,
            dump_every,
            [],
            os.path.join(path, "wall/wall_particles-"),
        )
    )

    sample_every = 2
    tdump_every = 10
    dump_every = int(tdump_every / dt)
    print("Velocity average - dump every :", dump_every)
    h = radius / 20
    bin_size = (domain[0], h, h)
    u.registerPlugins(
        mir.Plugins.createDumpAverage(
            "field",
            [pv_solvent],
            sample_every,
            dump_every,
            bin_size,
            ["velocities"],
            os.path.join(path, "solvent-vel-average"),
        )
    )

    # we can dump the wall sdf in xdmf + h5 format for visualization purpose
    # the dumped file is a grid with spacings h containing the SDF values
    u.dumpWalls2XDMF(
        [wall], h=(0.5, 0.5, 0.5), filename=os.path.join(path, "wall/wall-")
    )

    u.run(t_end, dt=dt)


def center_of_rbc(domain, radius, l, rc, n_rbc):
    """
    Génère une liste de centres pour placer des membranes (rbc) dans un cylindre.
    Les centres sont espacés régulièrement le long de l'axe x, et répartis en cercle dans la section yz.

    Args:
        domain (tuple): dimensions du domaine (l, L, L)
        radius (float): rayon du cylindre
        l (float): longueur du cylindre (axe x)
        rc (float): rayon de coupure (pour éviter les chevauchements)

    Returns:
        list: liste des centres sous forme de tuples (x, y, z)
    """
    list_center = []
    bound = 2
    center_yz = (domain[1] * 0.5, domain[2] * 0.5)
    print("Center_yz", center_yz)
    t = radius / 2  # rayon du cercle sur lequel placer les centres dans la section yz
    angles = np.array([0, 2 * np.pi / 3, 4 * np.pi / 3])  # 3 positions angulaires

    x = bound
    s = (l - 2 * bound) * 3 / n_rbc
    step = s

    while x < l - bound:
        y = center_yz[0] + t * np.sin(angles)
        z = center_yz[1] + t * np.cos(angles)
        for yi, zi in zip(y, z):
            list_center.append([x, yi, zi])
        x += step
    return list_center


def postprocess_cylinder_2D(
    path: str, h: float, R: float, Umax: float, path_cylinder=None
):
    import os
    import glob
    import h5py as h5
    import numpy as np
    import matplotlib.pyplot as plt

    files = sorted(glob.glob(os.path.join(path, "*.h5")))
    files = files[len(files) // 2 :]

    mean_v = None
    max_v_t = []
    for fname in files:
        with h5.File(fname, "r") as f:
            vel = np.array(f["velocities"])  
            print("Shape:", vel.shape)
            vel = vel.squeeze(axis=2)  
            v = vel[:, :, 0]  
            max_v_t.append(np.max(np.array(v)))
            if mean_v is None:
                mean_v = v
            else:
                mean_v += v
    mean_v /= len(files)

    mean_v_cy = None
    if path_cylinder is not None:
        files_cy = sorted(glob.glob(os.path.join(path_cylinder, "*.h5")))
        files_cy = files_cy[len(files_cy) // 2 :]
        max_v_t = []
        for fname in files_cy:
            with h5.File(fname, "r") as f:
                vel = np.array(f["velocities"])  
                print("Shape:", vel.shape)
                vel = vel.squeeze(axis=2)  
                v = vel[:, :, 0]  
                max_v_t.append(np.max(np.array(v)))
                if mean_v_cy is None:
                    mean_v_cy = v
                else:
                    mean_v_cy += v

        mean_v_cy /= len(files_cy)
    print("Shape mean_v: ", mean_v.shape)
    plt.figure(figsize=(12, 8))

    y = np.linspace(0, 2 * R, mean_v.shape[1])
    z = np.linspace(0, 2 * R, mean_v.shape[0])
    Y, Z = np.meshgrid(y, z)

    im = plt.contourf(Y, Z, mean_v, levels=50, cmap="RdYlBu_r")
    plt.colorbar(im, label="Vitesse u_x (m/s)")


    plt.xlabel("x (μm)")
    plt.ylabel("y (μm)")
    plt.title(f"Champ de vitesse u_x - Écoulement autour d'un cylindre (R={R}m)")
    plt.axis("equal")
    plt.grid(True, alpha=0.3)

    plt.figure(figsize=(12, 8))
    mean_v_norm = mean_v

    im2 = plt.contourf(Y, Z, mean_v_norm, levels=50, cmap="RdYlBu_r")
    plt.colorbar(im2, label="Vitesse u_x")

    circle = plt.Circle((R, R), R, color="black", fill=False)
    plt.gca().add_patch(circle)

    plt.xlabel("y (μm)")
    plt.ylabel("z (μm)")
    plt.axis("equal")
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(
        "fig/2D_cylinder_with_few_membranes_bis.png", dpi=100, bbox_inches="tight"
    )
    plt.close()
    radial_profile(mean_v, R, mean_v_cy)

    plt.close()

    plt.plot(max_v_t)
    plt.xlabel("simulation time")
    plt.ylabel("⟨v_z⟩_max")
    plt.title("Maximum of velocity")
    plt.legend()
    plt.grid()
    plt.savefig(
        "fig/vz_max_over_time_cylinder_with_few_membranes_bis_.png",
        dpi=100,
        bbox_inches="tight",
    )

    return mean_v


def radial_profile(mean_v, R, mean_v_cy=None):
    """
    Calcule le profil radial moyen de la vitesse
    mean_v: champ de vitesse 2D (45x45)
    R: rayon du cylindre
    """
    import numpy as np
    import matplotlib.pyplot as plt

    ny, nx = mean_v.shape

    x = np.linspace(0, 2 * R, nx)
    y = np.linspace(0, 2 * R, ny)

    x_center, y_center = R, R

    X, Y = np.meshgrid(x, y)
    r_grid = np.sqrt((X - x_center) ** 2 + (Y - y_center) ** 2)

    r_bins = np.linspace(0, R, 25)
    r_profile = []
    v_profile = []
    r_profile_cy = []
    v_profile_cy = []

    for i in range(len(r_bins) - 1):
        r1, r2 = r_bins[i], r_bins[i + 1]
        r_mid = (r1 + r2) / 2

        mask = (r_grid >= r1) & (r_grid < r2)

        if np.sum(mask) > 0:
            v_mean = np.mean(mean_v[mask])
            r_profile.append(r_mid)
            v_profile.append(v_mean)

    r_profile = np.array(r_profile)
    v_profile = np.array(v_profile)

    if mean_v_cy is not None:
        for i in range(len(r_bins) - 1):
            r1, r2 = r_bins[i], r_bins[i + 1]
            r_mid = (r1 + r2) / 2

            mask = (r_grid >= r1) & (r_grid < r2)

            if np.sum(mask) > 0:  
                v_mean_cy = np.mean(mean_v_cy[mask])
                r_profile_cy.append(r_mid)
                v_profile_cy.append(v_mean_cy)

        r_profile_cy = np.array(r_profile_cy)
        v_profile_cy = np.array(v_profile_cy)

   
    r_centered = r_profile

    v_max = np.max(v_profile)

    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    axs[0].plot(r_centered, v_profile, "+k", label="exp with RBCs")
    if mean_v_cy is not None:
        fit = np.polyfit(r_centered, v_profile_cy, deg=2)
        axs[0].plot(r_centered, v_profile_cy, "+b", label="exp without RBC")
        axs[0].plot(
            r_centered,
            np.polyval(fit, r_centered),
            "-k",
            label="fit on exp without RBC",
        )
    axs[0].axvline(x=R , color="r", linestyle="--", alpha=0.7, label="Cylinder surface")
    axs[0].set_xlabel("r (μm)")
    axs[0].set_ylabel("Mean velocity u_x (μm/s)")
    axs[0].set_title("Radial velocity")
    axs[0].grid(True, alpha=0.3)
    axs[0].legend()

    mean_v_slice_y = mean_v[:, mean_v.shape[1] // 2]

    y = np.linspace(0, 2 * R, mean_v.shape[0]) 
    z = np.linspace(0, 2 * R, mean_v.shape[1])

    axs[1].plot(y, mean_v_slice_y, "-o", color="green", label="exp with RBCs")
    if mean_v_cy is not None:
        mean_v_slice_y_cy = mean_v_cy[:, mean_v.shape[1] // 2]
        axs[1].plot(y, mean_v_slice_y_cy, "-o", color="blue", label="exp without RBC")
        fit = np.polyfit(y, mean_v_slice_y_cy, deg=2)
        axs[1].plot(y, np.polyval(fit, y), "-k", label="fit on exp without RBC")
    axs[1].set_xlabel("y (μm)")
    axs[1].set_ylabel("u_x (μm/s)")
    axs[1].set_title("Complete profile Along y axis")
    axs[1].grid(True, alpha=0.3)
    axs[1].legend()

    plt.tight_layout()
    plt.savefig(
        "fig/radial_profil_with_slice_cylinder_few_membranes_bis.png",
        dpi=100,
        bbox_inches="tight",
    )
    plt.close()

    return r_centered, v_profile


def parameters_cylinder(Vrbc, rc, Rrbc, Ht, l, R):

    # creation of the walls
    # we create a cylindrical pipe passing through the center of the domain along
    L = 2 * R + 2 * rc
    domain = (l, L, L)
    center = (domain[1] * 0.5, domain[2] * 0.5)  # center in the (yz) plane

    Vcy = np.pi * R**2 * l
    print("Volume du cylindre : ", Vcy)
    print("Volume RBC : ", Vrbc)
    n_rbc = Ht * Vcy / Vrbc
    print(f"Nombre de rbc pour {Ht} d'hematocrit :", n_rbc)
    h = R / 25
    list_center = center_of_rbc(domain, R, l, rc, n_rbc)
    nb_rbc_eff = len(list_center)
    print("Number of RBC :", len(list_center))
    print("Effective hematocrit :", nb_rbc_eff * Vrbc / Vcy)
    return list_center, domain, R, center


if __name__ == "__main__":
    Re = 0.1
    t_end = 6000000
    ranks = (1, 1, 1)
    verbose = True
    params = set_parameters(Re, verbose)
    path = "h5/cylinder_few_membranes_bis_Re_01_long/"
    radius = params["R"]
    Umax = params["Vmax"]
    h = radius / 20
    # cylinder_flow(ranks, path, params, t_end)
    path_cylinder = "h5/cylinder/compute_profile_L_h_h"
    mean_v = postprocess_cylinder_2D(path, h, radius, Umax)
