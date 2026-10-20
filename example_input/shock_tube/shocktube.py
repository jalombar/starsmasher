#!/usr/bin/env python3
"""A periodic shock tube in 3D, after Lombardi et al. (1999, J. Comput. Phys.
152, 687, sec. 5): two slabs of gamma-law gas at rest, side by side in a
periodic cell -1/2 < x <= 1/2 of cross-section 1 x 1, so that the cell holds
two interfaces and each launches a shock and a rarefaction.

  python3 shocktube.py [--case lowM|highM|gam3] [--n 10000] [--np 2] [--exe PATH]

Each slab is an hcp lattice with a random orientation, relaxed on its own
in a periodic half cell (nintvar=1 with a drag) to remove the defects of
the cut.  The relaxed slabs are then put side by side (make_ic.py combine)
and the shock tube runs to t=5 in the directory tube/.  Needs numpy and
scipy.  Every particle has the same mass, so the slabs differ in number
density.
"""
import argparse
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = {  # left state rho=1, p=1 in every case
    "lowM": dict(rho_r=0.25, p_r=5 / 2 ** (16 / 3), gam=5 / 3),
    "highM": dict(rho_r=0.25, p_r=1.2402e-3, gam=5 / 3),
    "gam3": dict(rho_r=0.25, p_r=8.78e-4, gam=3.0),
}


def write_inputs(d, **kw):
    lines = [f" {k}='{v}'," if isinstance(v, str) else
             f" {k}={v}," if isinstance(v, int) else f" {k}={float(v)!r}," for k, v in kw.items()]
    with open(os.path.join(d, "sph.input"), "w") as f:
        f.write(" &input\n" + "\n".join(lines).rstrip(",") + " /\n")
    with open(os.path.join(d, "sph.init"), "w") as f:
        f.write(" &initt\n iname='pbx' /\n")


def run(d, a):
    with open(os.path.join(d, "stdout.txt"), "w") as f:
        subprocess.run(["mpirun", "-np", str(a.np), a.exe], cwd=d, stdout=f,
                       stderr=subprocess.STDOUT, check=True)


def last_dump(d):
    return os.path.join(d, sorted(f for f in os.listdir(d)
                                  if f.startswith("out") and f.endswith(".sph"))[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="lowM", choices=CASES)
    ap.add_argument("--n", type=int, default=10000, help="particles in the cell")
    ap.add_argument("--np", type=int, default=2, help="MPI ranks")
    ap.add_argument("--exe", default=os.path.join(HERE, "..", "..", "parallel_bleeding_edge",
                                                   "parallel_bleeding_edge_cpu_sph"))
    a = ap.parse_args()
    a.exe = os.path.abspath(a.exe)
    c = CASES[a.case]
    g = c["gam"]
    common = dict(ngr=0, neos=0, gam=g, nnopt=39, startfile1="start.sph", boxy=1.0, boxz=1.0)
    nl = round(a.n / (1 + c["rho_r"]))
    slabs = {"l": (1.0, 1.0, nl), "r": (c["rho_r"], c["p_r"], a.n - nl)}
    for side, (rho, p, n) in slabs.items():
        d = os.path.join(HERE, a.case, "slab_" + side)
        os.makedirs(d, exist_ok=True)
        subprocess.run([sys.executable, os.path.join(HERE, "make_ic.py"), "slab",
                        os.path.join(d, "start.sph"), "--rho", str(rho), "--p", str(p),
                        "--gam", str(g), "--n", str(n), "--seed", str(1 if side == "l" else 2)],
                       check=True)
        # ts: sound crossing time of one particle spacing; drag for 240 of them
        ts = (np.sqrt(2) * 0.5 / n) ** (1 / 3) / np.sqrt(g * p / rho)
        write_inputs(d, **common, boxx=0.5, nintvar=1, nav=3, alpha=1.0, beta=2.0, nrelax=1,
                     trelax=3 * ts, treloff=240 * ts, tf=240 * ts, dtout=60 * ts)
        print(f"relaxing slab {side} ({n} particles) in {d}", flush=True)
        run(d, a)
    d = os.path.join(HERE, a.case, "tube")
    os.makedirs(d, exist_ok=True)
    subprocess.run([sys.executable, os.path.join(HERE, "make_ic.py"), "combine",
                    os.path.join(d, "start.sph"),
                    last_dump(os.path.join(HERE, a.case, "slab_l")),
                    last_dump(os.path.join(HERE, a.case, "slab_r")),
                    "--rho-l", "1", "--p-l", "1", "--rho-r", str(c["rho_r"]),
                    "--p-r", str(c["p_r"]), "--gam", str(g)], check=True)
    # alpha = 0.5 gam/2 and beta = gam/2 fit all three cases well
    write_inputs(d, **common, boxx=1.0, nintvar=2, nav=3, alpha=0.5 * g / 2, beta=g / 2,
                 nrelax=0, tf=5.0, dtout=0.1)
    print(f"running the shock tube in {d}", flush=True)
    run(d, a)
    print(f"done: dumps out*.sph in {d}, read them with sphio.read()")


if __name__ == "__main__":
    main()
