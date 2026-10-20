"""Read and write StarSmasher out*.sph dumps (ncooling=0 layout).

A dump is a Fortran sequential unformatted file: one header record, one
record per particle, and a trailing record holding ntot again.
"""
import numpy as np

HEADER = np.dtype([
    ("ntot", "<i4"), ("nnopt", "<i4"), ("hco", "<f8"), ("hfloor", "<f8"),
    ("sep0", "<f8"), ("tf", "<f8"), ("dtout", "<f8"), ("nout", "<i4"),
    ("nit", "<i4"), ("t", "<f8"), ("nav", "<i4"), ("alpha", "<f8"),
    ("beta", "<f8"), ("tjumpahead", "<f8"), ("ngr", "<i4"), ("nrelax", "<i4"),
    ("trelax", "<f8"), ("dt", "<f8"), ("omega2", "<f8"), ("ncooling", "<i4"),
    ("erad", "<f8"), ("ndisplace", "<i4"), ("displacex", "<f8"),
    ("displacey", "<f8"), ("displacez", "<f8"),
])
PART = np.dtype([
    ("x", "<f8"), ("y", "<f8"), ("z", "<f8"), ("am", "<f8"), ("hp", "<f8"),
    ("rho", "<f8"), ("vx", "<f8"), ("vy", "<f8"), ("vz", "<f8"),
    ("vxdot", "<f8"), ("vydot", "<f8"), ("vzdot", "<f8"), ("u", "<f8"),
    ("udot", "<f8"), ("grpot", "<f8"), ("meanmolecular", "<f8"),
    ("cc", "<i4"), ("divv", "<f8"),
])


def read(fn):
    raw = open(fn, "rb").read()
    pos = 0

    def record():
        nonlocal pos
        nb = int(np.frombuffer(raw, "<i4", 1, pos)[0])
        body = raw[pos + 4: pos + 4 + nb]
        pos += nb + 8
        return body

    hb = record()
    hdr = np.frombuffer(hb[:HEADER.itemsize], HEADER)[0]
    n = int(hdr["ntot"])
    p = np.empty(n, PART)
    reclen = PART.itemsize
    # the particle records are contiguous, each framed by 4-byte markers
    block = np.frombuffer(raw, np.uint8, n * (reclen + 8), pos)
    block = block.reshape(n, reclen + 8)[:, 4:4 + reclen].copy()
    p[:] = block.view(PART).ravel()
    pos += n * (reclen + 8)
    nchk = int(np.frombuffer(record(), "<i4", 1)[0])
    if nchk != n:
        raise ValueError(f"{fn}: trailing ntot {nchk} != {n}")
    return hdr, p


def write(fn, p, t=0.0, nnopt=23, alpha=1.0, beta=2.0):
    hdr = np.zeros(1, HEADER)
    hdr["ntot"] = len(p)
    hdr["nnopt"] = nnopt
    hdr["t"] = t
    hdr["nav"] = 3
    hdr["alpha"] = alpha
    hdr["beta"] = beta
    hdr["tjumpahead"] = 1e30
    hdr["trelax"] = 1e30
    with open(fn, "wb") as f:
        def rec(b):
            m = np.array([len(b)], "<i4").tobytes()
            f.write(m + b + m)
        rec(hdr.tobytes())
        pb = np.ascontiguousarray(p.astype(PART))
        for i in range(len(pb)):
            rec(pb[i:i + 1].tobytes())
        rec(np.array([len(p)], "<i4").tobytes())
