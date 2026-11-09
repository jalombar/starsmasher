#!/usr/bin/env python3
"""Read StarSmasher snapshots into numpy arrays: out*.sph or out*.h5.

Both formats give the same thing: a header (a dict) and a dict of numpy arrays,
one per column, named like the Fortran variables of the code.

    import StSm_snapshot as ss

    header, p = ss.read('out0003.h5')          # or 'out0003.sph'
    print(header['t'], p['rho'].max())

    # only some columns; from an HDF5 file the others are never read
    header, p = ss.read('out0003.h5', columns=['x', 'y', 'z', 'rho'])

The columns are x y z am hp rho vx vy vz vxdot vydot vzdot u udot grpot
meanmolecular cc divv, and ueq tthermal as well when the run had ncooling not 0.
Everything is float64 except cc, which is int32.  Everything is in code units
(see the page on units in the documentation); `code_units` gives the factors.

out*.sph is a Fortran sequential unformatted file written by gfortran, which is
little-endian on x86-64 and ARM, and the reader assumes that.  Reading it also
checks the markers of every record and the ntot that ends the file.  out*.h5
needs h5py, and only that format does.

From the command line the script prints the header and a summary of a file:

    python3 StSm_snapshot.py out0003.sph
"""
import os
import struct
import sys

import numpy as np

# The header of out*.sph, in the order dump (output.f) writes it, and the same
# names and types as the attributes of /header in out*.h5.
HEADER = np.dtype([
    ("ntot", "<i4"), ("nnopt", "<i4"), ("hco", "<f8"), ("hfloor", "<f8"),
    ("sep0", "<f8"), ("tf", "<f8"), ("dtout", "<f8"), ("nout", "<i4"),
    ("nit", "<i4"), ("t", "<f8"), ("nav", "<i4"), ("alpha", "<f8"),
    ("beta", "<f8"), ("tjumpahead", "<f8"), ("ngr", "<i4"), ("nrelax", "<i4"),
    ("trelax", "<f8"), ("dt", "<f8"), ("omega2", "<f8"), ("ncooling", "<i4"),
    ("erad", "<f8"), ("lastnumstars", "<i4"), ("epoteat", "<f8"),
    ("ekineat", "<f8"), ("einteat", "<f8"),
])

# The columns of a particle record, in file order.  ueq and tthermal follow
# divv only when ncooling is not 0.
COLUMNS = [
    ("x", "<f8"), ("y", "<f8"), ("z", "<f8"), ("am", "<f8"), ("hp", "<f8"),
    ("rho", "<f8"), ("vx", "<f8"), ("vy", "<f8"), ("vz", "<f8"),
    ("vxdot", "<f8"), ("vydot", "<f8"), ("vzdot", "<f8"), ("u", "<f8"),
    ("udot", "<f8"), ("grpot", "<f8"), ("meanmolecular", "<f8"),
    ("cc", "<i4"), ("divv", "<f8"),
]
COOLING_COLUMNS = [("ueq", "<f8"), ("tthermal", "<f8")]

HDF5_MAGIC = b"\x89HDF\r\n\x1a\n"


def columns_for(ncooling):
    """The columns a snapshot holds: two more when ncooling is not 0."""
    return COLUMNS + (COOLING_COLUMNS if ncooling != 0 else [])


def is_hdf5(path):
    with open(path, "rb") as f:
        return f.read(8) == HDF5_MAGIC


def _pick(columns, available):
    if columns is None:
        return list(available)
    columns = list(columns)
    missing = [c for c in columns if c not in available]
    if missing:
        raise KeyError("no such column: %s (the file has %s)"
                       % (", ".join(missing), ", ".join(available)))
    return columns


def read_native(path, columns=None):
    """Read out*.sph (or restartrad.sph): returns (header, data).

    Both layouts are handled: the one of ncooling=0 and the one with ueq and
    tthermal added, chosen from ncooling in the header and confirmed by the
    size of the file.  A record whose markers are wrong, a file of the wrong
    size or a final ntot that disagrees with the header raises ValueError.
    """
    with open(path, "rb") as f:
        raw = f.read(4)
        if len(raw) < 4 or struct.unpack("<i", raw)[0] != HEADER.itemsize:
            raise ValueError("%s: not a StarSmasher snapshot (the first record "
                             "should be %d bytes)" % (path, HEADER.itemsize))
        hdr = np.frombuffer(f.read(HEADER.itemsize), HEADER)[0]
        if struct.unpack("<i", f.read(4))[0] != HEADER.itemsize:
            raise ValueError("%s: the header record is not closed" % path)
        n = int(hdr["ntot"])
        ncool = int(hdr["ncooling"])
        size = os.path.getsize(path)
        start = 8 + HEADER.itemsize
        layout = None
        for nc in (ncool, 0, 1):          # the header first; then each layout
            cols = columns_for(nc)
            reclen = sum(np.dtype(t).itemsize for _, t in cols)
            if size == start + n * (reclen + 8) + 12:
                layout = cols
                break
        if layout is None:
            raise ValueError("%s: size %d does not fit ntot=%d with either "
                             "layout (ncooling=%d); truncated or corrupted?"
                             % (path, size, n, ncool))
        reclen = sum(np.dtype(t).itemsize for _, t in layout)
        rec = np.dtype([("m0", "<i4")] + layout + [("m1", "<i4")])
        f.seek(start)
        block = np.fromfile(f, dtype=rec, count=n)
        f.seek(start + n * (reclen + 8))
        tail = f.read(12)
    if len(block) != n:
        raise ValueError("%s: expected %d particle records, found %d"
                         % (path, n, len(block)))
    if not (np.all(block["m0"] == reclen) and np.all(block["m1"] == reclen)):
        raise ValueError("%s: a particle record has the wrong length marker"
                         % path)
    m0, nchk, m1 = struct.unpack("<3i", tail)
    if m0 != 4 or m1 != 4 or nchk != n:
        raise ValueError("%s: the final record holds ntot=%d, the header says "
                         "%d; the file was not written to the end"
                         % (path, nchk, n))
    header = {name: hdr[name].item() for name in HEADER.names}
    names = _pick(columns, [c for c, _ in layout])
    data = {c: np.ascontiguousarray(block[c]) for c in names}
    return header, data


def read_hdf5(path, columns=None):
    """Read out*.h5: returns (header, data), the header being all of /header.

    Besides the 25 entries of the header of out*.sph it holds format,
    format_version and the settings a reader needs (munit, runit, gravconst,
    tunit, nintvar, neos, ...).  Only the datasets in `columns` are read.
    """
    try:
        import h5py
    except ImportError as e:
        raise ImportError("reading out*.h5 needs h5py (%s)" % e)
    with h5py.File(path, "r") as f:
        if "header" not in f or "particles" not in f:
            raise ValueError("%s: no /header and /particles: not a "
                             "StarSmasher HDF5 snapshot" % path)
        header = {}
        for k, v in f["header"].attrs.items():
            header[k] = v.item() if hasattr(v, "item") else v
        grp = f["particles"]
        order = [c for c, _ in COLUMNS + COOLING_COLUMNS if c in grp]
        order += [c for c in grp if c not in order]
        names = _pick(columns, order)
        data = {c: grp[c][()] for c in names}
    return header, data


def read(path, columns=None):
    """Read a snapshot of either format: returns (header, data).

    `data` maps each column name to a numpy array of length ntot.  `columns`
    is a list of names to keep (all by default); for out*.h5 only those are
    read from the disk, while out*.sph has to be read whole.
    """
    if is_hdf5(path):
        return read_hdf5(path, columns)
    return read_native(path, columns)


def describe(path):
    """{column: (units, description)} from an out*.h5 file."""
    import h5py
    with h5py.File(path, "r") as f:
        return {c: (f["particles"][c].attrs["units"],
                    f["particles"][c].attrs["description"])
                for c in f["particles"]}


def code_units(header=None, munit=None, runit=None, gravconst=None):
    """What a code unit is worth in cgs, as a dict.

    G, the unit of mass and the unit of length are 1; the others follow from
    munit, runit and gravconst.  An HDF5 header carries them; for out*.sph, which
    does not, give them or get the defaults of the code (MESA solar mass and
    radius, CODATA 2018 G).
    """
    header = header or {}
    munit = munit or header.get("munit", 1.9884098706980504e33)
    runit = runit or header.get("runit", 6.957e10)
    gravconst = gravconst or header.get("gravconst", 6.67430e-08)
    tunit = (runit ** 3 / (gravconst * munit)) ** 0.5
    vunit = runit / tunit
    return dict(mass=munit, length=runit, time=tunit, velocity=vunit,
                acceleration=vunit / tunit, density=munit / runit ** 3,
                energy=gravconst * munit ** 2 / runit,
                specific_energy=vunit ** 2,
                pressure=gravconst * munit ** 2 / runit ** 4)


def to_dataframe(data, columns=None):
    """A pandas DataFrame with the names SPH tools expect.

    am becomes m (mass) and hp becomes h (smoothing length, the kernel reaching
    2h), so that sarracen.SarracenDataFrame(df) takes it directly together with
    x, y, z and rho.  The other columns keep their names.
    """
    import pandas as pd
    df = pd.DataFrame({k: v for k, v in data.items()
                       if columns is None or k in columns})
    return df.rename(columns={"am": "m", "hp": "h"})


def compare(path_a, path_b):
    """Compare two snapshots bit for bit, whatever their formats.

    Returns (ndataset, datasets_differing, attributes_differing): how many
    columns both files hold, the names of those that differ in any bit (NaN
    included), and the names of the entries of the common 25-entry header that
    differ.  Floats are compared through their bit patterns.
    """
    ha, da = read(path_a)
    hb, db = read(path_b)
    common = [c for c in da if c in db]
    bad = [c for c in common if not _same_bits(da[c], db[c])]
    for c in set(da) ^ set(db):
        bad.append(c)               # a column only one of them holds
    badh = [k for k in HEADER.names
            if not _same_bits(np.asarray(ha[k], HEADER[k]),
                              np.asarray(hb[k], HEADER[k]))]
    return len(common), sorted(bad), badh


def _same_bits(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.dtype != b.dtype or a.shape != b.shape:
        return False
    if a.dtype.kind == "f":
        a, b = a.view("<u%d" % a.dtype.itemsize), b.view("<u%d" % b.dtype.itemsize)
    return bool(np.array_equal(a, b))


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    header, p = read(argv[1])
    fmt = "HDF5" if is_hdf5(argv[1]) else "native"
    print("%s (%s): %d particles, t = %.8g, nout = %d, nit = %d"
          % (argv[1], fmt, header["ntot"], header["t"], header["nout"],
             header["nit"]))
    for k in HEADER.names:
        print("  %-13s %s" % (k, header[k]))
    print("columns:")
    for c, a in p.items():
        print("  %-14s %-8s min % .6e  max % .6e" % (c, a.dtype, a.min(), a.max()))
    print("  total mass     % .12e" % p["am"].sum())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
