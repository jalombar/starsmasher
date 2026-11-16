#!/usr/bin/env python3
"""Two plots from one StarSmasher snapshot, out*.h5 or out*.sph.

    python3 StSm_snapshot_plots.py out0003.h5 plots/star

writes plots/star_rho_profile.png and plots/star_column_density.png:

  * the density of every particle against its distance from the centre of
    mass, with the median in logarithmic bins on top;
  * the column density on the x-y plane, as a two-dimensional histogram of the
    particle masses (the mass of the particles in a pixel over its area).

The units are the code units unless --cgs is given, which uses munit and runit
from the header of an HDF5 file (the defaults of the code for a native one, or
those of --munit and --runit).  This is an example to copy and change; the
reading is done by StSm_snapshot.py, which has to be in the same directory.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import StSm_snapshot as ss


def centre_of_mass(p):
    m = p["am"]
    return np.array([np.sum(m * p[k]) / np.sum(m) for k in ("x", "y", "z")])


def binned_median(r, y, nbins=40):
    edges = np.logspace(np.log10(r[r > 0].min()), np.log10(r.max()), nbins + 1)
    idx = np.digitize(r, edges) - 1
    centre, med = [], []
    for i in range(nbins):
        sel = y[idx == i]
        if len(sel) >= 5:
            centre.append(np.sqrt(edges[i] * edges[i + 1]))
            med.append(np.median(sel))
    return np.array(centre), np.array(med)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("snapshot", help="out*.h5 or out*.sph")
    ap.add_argument("prefix", help="the PNG files are <prefix>_*.png")
    ap.add_argument("--cgs", action="store_true", help="plot in g/cm^3, cm and g/cm^2")
    ap.add_argument("--munit", type=float, help="grams in the unit of mass (native files)")
    ap.add_argument("--runit", type=float, help="cm in the unit of length (native files)")
    ap.add_argument("--bins", type=int, help="pixels per side of the map (default: about 60 particles "
                    "per pixel, between 48 and 400)")
    ap.add_argument("--width", type=float, help="full width of the map, in the unit of "
                    "length (default: holds 99.5%% of the particles)")
    a = ap.parse_args()

    header, p = ss.read(a.snapshot, columns=["x", "y", "z", "am", "rho"])
    cu = ss.code_units(header, munit=a.munit, runit=a.runit)
    lscale, rscale = (cu["length"], cu["density"]) if a.cgs else (1.0, 1.0)
    sscale = cu["mass"] / cu["length"] ** 2 if a.cgs else 1.0
    lunit = "cm" if a.cgs else "runit"
    title = "%s   t = %.5g   N = %d" % (os.path.basename(a.snapshot), header["t"],
                                        header["ntot"])
    os.makedirs(os.path.dirname(os.path.abspath(a.prefix)), exist_ok=True)

    com = centre_of_mass(p)
    dx, dy, dz = p["x"] - com[0], p["y"] - com[1], p["z"] - com[2]
    r = np.sqrt(dx * dx + dy * dy + dz * dz)

    # density profile
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    keep = (r > 0) & (p["rho"] > 0)
    step = max(1, keep.sum() // 50000)            # at most ~50000 points on the page
    ax.scatter(r[keep][::step] * lscale, p["rho"][keep][::step] * rscale,
               s=1, alpha=0.25, color="tab:blue", rasterized=True,
               label="particles")
    rc, med = binned_median(r[keep], p["rho"][keep])
    ax.plot(rc * lscale, med * rscale, color="k", lw=1.6, label="median per bin")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("r from the centre of mass [%s]" % lunit)
    ax.set_ylabel("density [%s]" % ("g cm$^{-3}$" if a.cgs else "munit / runit$^3$"))
    ax.set_title(title, fontsize=9)
    ax.legend()
    fig.tight_layout()
    f1 = a.prefix + "_rho_profile.png"
    fig.savefig(f1, dpi=130)
    plt.close(fig)

    # column density on x-y, projected along z
    half = 0.5 * a.width if a.width else 1.05 * np.quantile(
        np.maximum(np.abs(dx), np.abs(dy)), 0.995)
    rng = [[-half, half], [-half, half]]
    # a histogram needs several particles per pixel, or it shows the lattice
    # the particles were laid out on instead of the star
    bins = a.bins or int(min(400, max(48, np.sqrt(len(r) / 60.0))))
    mass, _, _ = np.histogram2d(dx, dy, bins=bins, range=rng, weights=p["am"])
    area = (2 * half / bins) ** 2
    sigma = mass / area
    fig, ax = plt.subplots(figsize=(6.0, 5.2))
    with np.errstate(divide="ignore"):
        img = np.log10(np.where(sigma > 0, sigma * sscale, np.nan))
    im = ax.imshow(img.T, origin="lower", cmap="inferno",
                   extent=[-half * lscale, half * lscale, -half * lscale, half * lscale])
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("log10 column density [%s]" % ("g cm$^{-2}$" if a.cgs
                                               else "munit / runit$^2$"))
    ax.set_xlabel("x - x$_{cm}$ [%s]" % lunit)
    ax.set_ylabel("y - y$_{cm}$ [%s]" % lunit)
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    f2 = a.prefix + "_column_density.png"
    fig.savefig(f2, dpi=130)
    plt.close(fig)
    print("wrote", f1)
    print("wrote", f2)


if __name__ == "__main__":
    main()
