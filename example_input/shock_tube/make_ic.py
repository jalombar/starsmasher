"""Initial conditions for the periodic shock tube of Lombardi et al. (1999,
sec. 5): a cell -1/2 < x <= 1/2 (eq. 29) of cross-section 1 x 1, periodic
on all six faces, left slab on -1/2 < x <= 0 and right slab on 0 < x <= 1/2.
All particles have the same mass, so the slabs differ in number density.

  make_ic.py slab OUT --rho R --p P --gam G --n N [--seed S]
      One slab, 1/2 x 1 x 1, centred on the origin: an hcp lattice with a
      random orientation, cut to the box with its spacing tuned so that
      exactly N lattice points fall inside.  The cut leaves defects at the
      periodic faces, which the relaxation (nintvar=1, drag) removes.
      u holds the entropic variable A=p/rho^gam.

  make_ic.py combine OUT LEFT RIGHT --rho-l --p-l --rho-r --p-r --gam
      Put a relaxed left slab on -1/2 < x < 0 and a relaxed right slab on
      0 < x < 1/2.  All masses are scaled by one common factor so that the
      mean SPH density of the left slab is rho_l (the kernel estimate runs
      slightly high), velocities are zeroed, and u = A rho^(gam-1)/(gam-1)
      from each particle's own relaxed density, for nintvar=2.
"""
import argparse
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.spatial.transform import Rotation
import sphio

SLAB = np.array([0.5, 1.0, 1.0])          # one slab
LX, LY, LZ = 1.0, 1.0, 1.0                 # the whole periodic cell


def hcp_block(a, half):
    """hcp points of spacing a filling the cube |x|,|y|,|z| < half."""
    dy, dz = a * np.sqrt(3) / 2, a * np.sqrt(2 / 3)
    ni, nj, nk = (int(np.ceil(half / s)) + 2 for s in (a, dy, dz))
    i, j, k = np.meshgrid(np.arange(-ni, ni + 1), np.arange(-nj, nj + 1),
                          np.arange(-nk, nk + 1), indexing="ij")
    i, j, k = i.ravel(), j.ravel(), k.ravel()
    x = (i + 0.5 * (j % 2) + 0.5 * (k % 2)) * a
    y = (j + (k % 2) / 3.0) * dy
    z = k * dz
    return np.column_stack([x, y, z])


def rotated_lattice(n, rng):
    """exactly n points of a randomly oriented hcp lattice inside SLAB."""
    rot = Rotation.random(random_state=rng)
    shift = rng.uniform(-0.5, 0.5, 3)
    a0 = (np.sqrt(2) * np.prod(SLAB) / n) ** (1 / 3)   # hcp: n = sqrt(2)/a^3
    half = 0.5 * np.linalg.norm(SLAB) + 2 * a0
    def inside(a):
        p = rot.apply(hcp_block(a, half) + shift * a)
        return p[np.all(np.abs(p) < SLAB / 2, axis=1)]
    lo, hi = 0.97 * a0, 1.03 * a0
    for _ in range(60):
        a = 0.5 * (lo + hi)
        p = inside(a)
        if len(p) == n:
            return p, a
        lo, hi = (a, hi) if len(p) > n else (lo, a)
    # the count can step past n; drop the surplus at random
    p = inside(lo)
    return p[np.sort(rng.choice(len(p), n, replace=False))], lo


def slab(a):
    rng = np.random.default_rng(a.seed)
    xyz, spacing = rotated_lattice(a.n, rng)
    p = np.zeros(a.n, sphio.PART)
    p["x"], p["y"], p["z"] = xyz.T
    p["am"] = a.rho * np.prod(SLAB) / a.n
    # initial guess for h; StarSmasher solves for it, but its periodic
    # padding (3 h_max) must already cover the first neighbour search
    p["hp"] = spacing * a.hp_factor
    p["rho"] = a.rho
    p["u"] = a.p / a.rho ** a.gam
    p["meanmolecular"] = 1.0
    sphio.write(a.out, p)
    print(f"{a.out}: n={a.n} m={p['am'][0]:.4e} spacing={spacing:.5f} A={p['u'][0]:.5g}")


def combine(a):
    _, L = sphio.read(a.left)
    _, R = sphio.read(a.right)
    L = L.copy(); R = R.copy()
    L["x"] -= 0.25
    R["x"] += 0.25
    f = a.rho_l / L["rho"].mean()
    for q, rho, pr, cc in ((L, a.rho_l, a.p_l, 0), (R, a.rho_r, a.p_r, 1)):
        q["am"] *= f
        q["rho"] *= f
        q["u"] = pr / rho ** a.gam * q["rho"] ** (a.gam - 1) / (a.gam - 1)
        q["cc"] = cc
        print(f"  slab {cc}: mean rho {q['rho'].mean():.5f} (nominal {rho:g}), "
              f"scatter {q['rho'].std() / q['rho'].mean():.2e}")
    p = np.concatenate([L, R])
    for k in ("vx", "vy", "vz", "vxdot", "vydot", "vzdot", "udot", "divv"):
        p[k] = 0
    sphio.write(a.out, p)
    print(f"{a.out}: n={len(p)}, mass factor {f:.5f}, total mass {p['am'].sum():.5f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("slab")
    s.add_argument("out")
    for k in ("rho", "p", "gam"):
        s.add_argument("--" + k, type=float, required=True)
    s.add_argument("--n", type=int, required=True)
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--hp-factor", type=float, default=1.0,
                   help="initial h in units of the lattice spacing (1 suits nnopt=39)")
    c = sub.add_parser("combine")
    c.add_argument("out"); c.add_argument("left"); c.add_argument("right")
    for k in ("rho-l", "p-l", "rho-r", "p-r", "gam"):
        c.add_argument("--" + k, type=float, required=True)
    a = ap.parse_args()
    slab(a) if a.cmd == "slab" else combine(a)
