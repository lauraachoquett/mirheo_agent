#!/usr/bin/env python

import mirheo as mir
from dpdprops import create_dpd_params_from_Re_Ma
import numpy as np
import os
import h5py as h5
import glob
import matplotlib.pyplot as plt
import json
from datetime import datetime


def cylinder_flow(
    Umax, rc, nd, ranks, L, domain, mass, Re, Ma, center, radius, h, path, t_end
):
    params = create_dpd_params_from_Re_Ma(
        Re=Re, Ma=Ma, U=Umax, L=2 * radius, nd=nd, rc=rc, mass=mass
    )

    dt = params.get_max_dt()
    nu = params.kinematic_viscosity()
    mu = params.dynamic_viscosity()
    print(f" Nu : {nu}, Mu : {mu}, Re : {Re}, Ma : {Ma} and dt :{dt}")

    f = Umax * 4 * nu * mass / radius**2
    print("F :", f)
    time = datetime.now().isoformat()
    params_save = {
        "Umax": Umax,
        "nd": nd,
        "rc": rc,
        "dt": dt,
        "ranks": ranks,
        "L": L,
        "domain": domain,
        "mass": mass,
        "Re": Re,
        "Ma": Ma,
        "t_end": t_end,
        "center": center,
        "radius": radius,
        "h": h,
        "force": f,
        "nu": nu,
        "timestamp": time,
    }
    os.makedirs("params", exist_ok=True)
    with open(f"params/params_cylinder_{time}.json", "w") as file:
        json.dump(params_save, file, indent=4)

    u = mir.Mirheo(ranks, domain, debug_level=0, log_filename="log")

    pv = mir.ParticleVectors.ParticleVector("pv", mass=1.0)
    ic = mir.InitialConditions.Uniform(nd)
    dpd = mir.Interactions.Pairwise("dpd", rc, kind="DPD", **params.to_interactions())
    vv = mir.Integrators.VelocityVerlet_withConstForce("vv", force=(f, 0, 0))

    u.registerInteraction(dpd)
    u.registerParticleVector(pv, ic)
    u.registerIntegrator(vv)

    wall = mir.Walls.Cylinder(
        "cylinder", center=center, radius=radius, axis="x", inside=True
    )

    u.registerWall(wall)  # register the wall in the coordinator

    # we now create the frozen particles of the walls
    # the following command is running a simulation of a solvent with given density equilibrating with dpd interactions and vv integrator
    # for 1000 steps.
    # It then selects the frozen particles according to the walls geometry, register and returns the newly created particle vector.
    pv_frozen = u.makeFrozenWallParticles(
        pvName="wall",
        walls=[wall],
        interactions=[dpd],
        integrator=vv,
        number_density=nd,
        dt=dt,
    )

    # set the wall for pv
    # this is required for non-penetrability of the solvent thanks to bounce-back
    # this will also remove the initial particles which are not inside the wall geometry
    u.setWall(wall, pv)

    # now the pv also interacts with the frozen particles!
    u.setInteraction(dpd, pv, pv)
    u.setInteraction(dpd, pv, pv_frozen)

    # pv_frozen do not move, only pv needs an integrator in this case
    u.setIntegrator(vv, pv)

    u.registerPlugins(mir.Plugins.createStats("stats", every=5000))
    dump_every = 5000
    u.registerPlugins(
        mir.Plugins.createDumpParticles(
            "part_dump", pv, dump_every, [], "h5/cylinder_bis/solvent_particles-"
        )
    )

    # we can also dump the frozen particles for visualization purpose
    u.registerPlugins(
        mir.Plugins.createDumpParticles(
            "part_dump_wall",
            pv_frozen,
            dump_every,
            [],
            "h5/cylinder_bis/wall_particles-",
        )
    )

    sample_every = 2
    tdump_every = 10
    dump_every = int(tdump_every / dt)
    h = radius / 20
    bin_size = (domain[0], h, h)
    u.registerPlugins(
        mir.Plugins.createDumpAverage(
            "field",
            [pv],
            sample_every,
            dump_every,
            bin_size,
            ["velocities"],
            os.path.join(path, "solvent-"),
        )
    )

    # we can dump the wall sdf in xdmf + h5 format for visualization purpose
    # the dumped file is a grid with spacings h containing the SDF values
    u.dumpWalls2XDMF([wall], h=(0.5, 0.5, 0.5), filename="h5/cylinder/wall")

    u.run(t_end, dt=dt)


def postprocess_cylinder_2D(path: str, h: float, R: float, Umax: float):
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
            vel = np.array(f["velocities"])  # shape: (45, 45, 1, 3)
            vel = vel.squeeze(axis=2)  # shape: (45, 45, 3)
            v = vel[:, :, 0]  # vitesse selon x
            max_v_t.append(np.max(np.array(v)))
            if mean_v is None:
                mean_v = v
            else:
                mean_v += v

    mean_v /= len(files)
    print("Shape mean_v: ", mean_v.shape)
    # Création du plot 2D
    plt.figure(figsize=(12, 8))

    # Coordonnées pour le plot (assumant un domaine centré)
    y = np.linspace(0, 2 * R, mean_v.shape[1])
    z = np.linspace(0, 2 * R, mean_v.shape[0])
    Y, Z = np.meshgrid(y, z)

    # Plot de la vitesse avec colormap
    im = plt.contourf(Y, Z, mean_v, levels=50, cmap="RdYlBu_r")
    plt.colorbar(im, label="Vitesse u_x (m/s)")

    # Lignes de courant pour visualiser l'écoulement
    plt.streamplot(
        Y, Z, mean_v, vel[:, :, 1], density=1.5, color="white", linewidth=0.8
    )

    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.title(f"Champ de vitesse u_x - Écoulement autour d'un cylindre (R={R}m)")
    plt.axis("equal")
    plt.grid(True, alpha=0.3)

    # Normalisation par Umax pour adimensionner
    plt.figure(figsize=(12, 8))
    mean_v_norm = mean_v

    im2 = plt.contourf(Y, Z, mean_v_norm, levels=50, cmap="RdYlBu_r")
    plt.colorbar(im2, label="Vitesse u_x")

    # Ajout du cylindre (cercle de rayon R centré à l'origine)
    circle = plt.Circle((R, R), R, color="black", fill=False)
    plt.gca().add_patch(circle)

    plt.xlabel("y (m)")
    plt.ylabel("z (m)")
    plt.title(f"Champ de vitesse adimensionnel ")
    plt.axis("equal")
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("fig/2D_cylinder.png", dpi=100, bbox_inches="tight")
    plt.close()
    radial_profile(mean_v, R, Umax)

    plt.close()

    plt.plot(max_v_t)
    plt.xlabel("simulation time")
    plt.ylabel("⟨v_z⟩_max")
    plt.title("Maximum of velocity")
    plt.legend()
    plt.grid()
    plt.savefig("fig/vz_max_over_time_cylinder.png", dpi=100, bbox_inches="tight")

    return mean_v


def radial_profile(mean_v, R, Umax):
    """
    Calcule le profil radial moyen de la vitesse
    mean_v: champ de vitesse 2D (45x45)
    R: rayon du cylindre
    """
    import numpy as np

    import matplotlib.pyplot as plt

    ny, nx = mean_v.shape

    # Coordonnées du domaine [0, 2R] x [0, 2R]
    x = np.linspace(0, 2 * R, nx)
    y = np.linspace(0, 2 * R, ny)

    # Centre du cylindre
    x_center, y_center = R, R

    # Calcul des distances radiales depuis le centre
    X, Y = np.meshgrid(x, y)
    r_grid = np.sqrt((X - x_center) ** 2 + (Y - y_center) ** 2)

    # Définition des rayons pour le profil (de 0 à R avec pas adaptatif)
    r_bins = np.linspace(0, R, 25)  # 25 points radiaux
    r_profile = []
    v_profile = []

    for i in range(len(r_bins) - 1):
        r1, r2 = r_bins[i], r_bins[i + 1]
        r_mid = (r1 + r2) / 2

        # Masque pour la couronne radiale
        mask = (r_grid >= r1) & (r_grid < r2)

        if np.sum(mask) > 0:  # Si des points existent dans cette couronne
            v_mean = np.mean(mean_v[mask])
            r_profile.append(r_mid)
            v_profile.append(v_mean)

    r_profile = np.array(r_profile)
    v_profile = np.array(v_profile)

    # Transposition pour avoir r ∈ [-R, R]
    r_centered = r_profile
    fit = np.polyfit(r_centered, v_profile, deg=2)

    v_max = np.max(v_profile)

    # Ajout d'un deuxième graphique avec subplot pour mean_v_slice
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    # Profil radial
    axs[0].plot(r_centered, v_profile, "+k", label="exp")
    axs[0].plot(r_centered, Umax * (1 - r_centered**2 / R**2), "-b", label="théorie")
    axs[0].plot(r_centered, np.polyval(fit, r_centered), "-k", label="fit")
    axs[0].axvline(x=R, color="r", linestyle="--", alpha=0.7, label="Surface cylindre")
    axs[0].set_xlabel("r (m)")
    axs[0].set_ylabel("Mean velocity u_x (m/s)")
    axs[0].set_title("Radial velocity")
    axs[0].grid(True, alpha=0.3)
    axs[0].legend()

    mean_v_slice_y = mean_v[:, mean_v.shape[1] // 2]
    y_center = np.linspace(-R, R, mean_v.shape[0])
    z = np.linspace(0, 2 * R, mean_v.shape[1])
    fit = np.polyfit(y_center, mean_v_slice_y, deg=2)

    axs[1].plot(y_center, mean_v_slice_y, "+k", label="exp")
    axs[1].plot(y_center, Umax * (1 - y_center**2 / R**2), "-b", label="théorie")

    axs[1].plot(y_center, np.polyval(fit, y_center), "-k", label="fit")
    axs[1].set_xlabel("y (m)")
    axs[1].set_ylabel("u_x (m/s)")
    axs[1].set_title("Complete profile Along y axis")
    axs[1].grid(True, alpha=0.3)
    axs[1].legend()

    plt.tight_layout()
    plt.savefig("fig/radial_profil_with_slice.png", dpi=100, bbox_inches="tight")
    plt.close()

    return r_centered, v_profile


if __name__ == "__main__":
    Umax = 1
    nd = 10.0
    rc = (10 / nd) ** (1 / 3)
    ranks = (1, 1, 1)

    L = 16
    diameter_rbc = 5
    Vrbc = 31.157

    l = 12 * diameter_rbc
    L = 16
    domain = (l, L, L)
    mass = 1
    Re = 0.1
    Ma = 0.02
    t_end = 1000000

    # Save parameters to JSON with timestamp

    path = "h5/cylinder/compute_profile_L_h_h"
    # creation of the walls
    # we create a cylindrical pipe passing through the center of the domain along
    center = (domain[1] * 0.5, domain[2] * 0.5)  # center in the (yz) plane
    radius = 0.5 * domain[1] - rc  # radius needs to be smaller than half of the domain
    # because of the frozen particles
    h = radius / 25
    print("R :", radius)
    print("2R/h :", 2 * radius / h)

    # cylinder_flow(
    #     Umax, rc, nd, ranks, L, domain, mass, Re, Ma, center, radius, h, path, t_end
    # )
    mean_v = postprocess_cylinder_2D(path, h, radius, Umax)
