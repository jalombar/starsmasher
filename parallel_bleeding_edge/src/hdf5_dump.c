/*
 * hdf5_dump.c -- writes a StarSmasher snapshot as an HDF5 file.
 *
 * Built only with `make HDF5=1`; without it hdf5_stub.f takes its place and
 * says the writer is not there.  dumph5 in output.f calls these functions,
 * through iso_c_binding, from the same place and with the same arrays that
 * dump uses for out*.sph, so the two files hold the same numbers.  Only rank 0
 * writes.  The C API is used, and not the Fortran modules of HDF5, because the
 * .mod files depend on the version of the Fortran compiler and the C API does
 * not.
 *
 * File layout, version 1 (docs/source/using/hdf5.rst has the full description):
 *
 *   /header      attributes: the 25 entries of the header of out*.sph, with
 *                the same names and types (int32 for integer, float64 for
 *                real*8); then format, format_version and the settings of
 *                sph.input that are needed to read the file (units, nintvar,
 *                neos, ...)
 *   /particles   one dataset of length ntot for every column of out*.sph,
 *                named like the Fortran variable (float64; cc is int32), each
 *                with the attributes `description` and `units`
 *
 * The datasets are not compressed: measured on real data, deflate (levels 1
 * and 4, with the shuffle filter) saves about a quarter of the size and makes
 * the write about eleven times slower and the read eight to nine times slower.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

#define FORMAT_VERSION 1

typedef struct {
    const char *name;
    int isint;
} entry;

/* The header of out*.sph, in the order of the Fortran write in dump. */
#define NHDR 25
static const entry HDR[NHDR] = {
    {"ntot", 1}, {"nnopt", 1}, {"hco", 0}, {"hfloor", 0}, {"sep0", 0},
    {"tf", 0}, {"dtout", 0}, {"nout", 1}, {"nit", 1}, {"t", 0}, {"nav", 1},
    {"alpha", 0}, {"beta", 0}, {"tjumpahead", 0}, {"ngr", 1}, {"nrelax", 1},
    {"trelax", 0}, {"dt", 0}, {"omega2", 0}, {"ncooling", 1}, {"erad", 0},
    {"lastnumstars", 1}, {"epoteat", 0}, {"ekineat", 0}, {"einteat", 0}
};

/* Settings of sph.input needed to read the file, in the order dumph5 fills
 * them.  tunit is derived: sqrt(runit^3/(gravconst*munit)). */
#define NEXT 17
static const entry EXT[NEXT] = {
    {"munit", 0}, {"runit", 0}, {"gravconst", 0}, {"tunit", 0},
    {"nintvar", 1}, {"neos", 1}, {"nkernel", 1}, {"gflag", 1},
    {"nusegpus", 1}, {"nselfgravity", 1}, {"nblock", 1}, {"dynhco", 1},
    {"gam", 0}, {"reat", 0}, {"boxx", 0}, {"boxy", 0},
    {"boxz", 0}
};

typedef struct {
    const char *name;
    const char *description;
    const char *units;
} column;

/* Units are in terms of munit, runit and tunit (see /header), with G = 1 in
 * code units.  The internal energy variable u is a specific energy only for
 * nintvar=2; /header/nintvar says which it is. */
static const column COLS[] = {
    {"x", "x coordinate of the particle", "runit"},
    {"y", "y coordinate of the particle", "runit"},
    {"z", "z coordinate of the particle", "runit"},
    {"am", "mass of the particle", "munit"},
    {"hp", "smoothing length h (the kernel reaches 2h)", "runit"},
    {"rho", "density", "munit/runit^3"},
    {"vx", "x velocity", "runit/tunit"},
    {"vy", "y velocity", "runit/tunit"},
    {"vz", "z velocity", "runit/tunit"},
    {"vxdot", "x acceleration", "runit/tunit^2"},
    {"vydot", "y acceleration", "runit/tunit^2"},
    {"vzdot", "z acceleration", "runit/tunit^2"},
    {"u", "internal energy variable: specific internal energy for nintvar=2, "
          "entropic variable A for nintvar=1, ln A for nintvar=3",
     "(runit/tunit)^2 for nintvar=2"},
    {"udot", "time derivative of u", "(runit/tunit)^2/tunit for nintvar=2"},
    {"grpot", "gravitational potential at the particle (energy per unit mass)",
     "(runit/tunit)^2"},
    {"meanmolecular", "mean molecular weight times the atomic mass unit",
     "g"},
    {"cc", "integer label of the component (star) the setup routines gave "
           "the particle", "dimensionless"},
    {"divv", "divergence of the velocity", "1/tunit"},
    {"ueq", "equilibrium specific internal energy the cooling relaxes "
            "towards (written when ncooling is not 0)", "(runit/tunit)^2"},
    {"tthermal", "thermal timescale of the cooling (written when ncooling "
                 "is not 0)", "tunit"},
};
#define NCOLS ((int)(sizeof(COLS) / sizeof(COLS[0])))

static struct {
    hid_t file, header, particles;
    hsize_t ntot;
    char path[1024], tmp[1032];
} S = {-1, -1, -1, 0, "", ""};

int ss_h5_available(void) { return 1; }

static int fail(const char *what)
{
    fprintf(stderr, " hdf5_dump: %s failed\n", what);
    if (S.particles >= 0) H5Gclose(S.particles);
    if (S.header >= 0) H5Gclose(S.header);
    if (S.file >= 0) H5Fclose(S.file);
    S.particles = S.header = S.file = -1;
    if (S.tmp[0]) remove(S.tmp);
    return 1;
}

static int attr_scalar(hid_t loc, const char *name, hid_t ftype,
                       hid_t mtype, const void *value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t a;
    herr_t r;
    if (space < 0) return 1;
    a = H5Acreate2(loc, name, ftype, space, H5P_DEFAULT, H5P_DEFAULT);
    if (a < 0) { H5Sclose(space); return 1; }
    r = H5Awrite(a, mtype, value);
    H5Aclose(a);
    H5Sclose(space);
    return r < 0;
}

static int attr_string(hid_t loc, const char *name, const char *text)
{
    hid_t t = H5Tcopy(H5T_C_S1);
    int bad;
    if (t < 0) return 1;
    H5Tset_size(t, H5T_VARIABLE);
    H5Tset_cset(t, H5T_CSET_UTF8);
    bad = attr_scalar(loc, name, t, t, &text);
    H5Tclose(t);
    return bad;
}

static int attr_int(hid_t loc, const char *name, int v)
{
    return attr_scalar(loc, name, H5T_STD_I32LE, H5T_NATIVE_INT, &v);
}

static int attr_double(hid_t loc, const char *name, double v)
{
    return attr_scalar(loc, name, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, &v);
}

static int attr_entry(hid_t loc, const entry *e, double v)
{
    /* integers travel as doubles, exactly, since they fit in 53 bits */
    if (e->isint) return attr_int(loc, e->name, (int)v);
    return attr_double(loc, e->name, v);
}

/* Opens <fname>.tmp and writes /header.  hdr has the 25 entries of the header
 * of out*.sph in file order, ext the NEXT settings; integers are passed as
 * doubles.  The file gets its final name in ss_h5_close, so that a run killed
 * in the middle of a write never leaves a half-written out*.h5 behind. */
int ss_h5_create(const char *fname, int ntot, const double *hdr, int nhdr,
                 const double *ext, int next)
{
    int i;
    unsigned maj, min, rel;
    char ver[64];
    hid_t fapl;

    if (nhdr != NHDR || next != NEXT) {
        fprintf(stderr, " hdf5_dump: header has %d/%d entries, expected "
                "%d/%d\n", nhdr, next, NHDR, NEXT);
        return 1;
    }
    /* One writer, one file: the locking of HDF5 only gets in the way, and
     * fails outright on some network file systems. */
    setenv("HDF5_USE_FILE_LOCKING", "FALSE", 0);

    S.ntot = (hsize_t)ntot;

    if (strlen(fname) + 5 > sizeof(S.tmp)) return fail("file name too long");
    strcpy(S.path, fname);
    snprintf(S.tmp, sizeof(S.tmp), "%s.tmp", fname);

    fapl = H5P_DEFAULT;
    S.file = H5Fcreate(S.tmp, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    if (S.file < 0) return fail("H5Fcreate");
    S.header = H5Gcreate2(S.file, "/header", H5P_DEFAULT, H5P_DEFAULT,
                          H5P_DEFAULT);
    if (S.header < 0) return fail("creating /header");
    S.particles = H5Gcreate2(S.file, "/particles", H5P_DEFAULT, H5P_DEFAULT,
                             H5P_DEFAULT);
    if (S.particles < 0) return fail("creating /particles");

    for (i = 0; i < NHDR; i++)
        if (attr_entry(S.header, &HDR[i], hdr[i]))
            return fail("writing a header attribute");
    if (attr_string(S.header, "format", "StarSmasher HDF5") ||
        attr_int(S.header, "format_version", FORMAT_VERSION))
        return fail("writing format attributes");
    for (i = 0; i < NEXT; i++)
        if (attr_entry(S.header, &EXT[i], ext[i]))
            return fail("writing a settings attribute");
    H5get_libversion(&maj, &min, &rel);
    snprintf(ver, sizeof(ver), "%u.%u.%u", maj, min, rel);
    if (attr_string(S.header, "hdf5_library_version", ver))
        return fail("writing the library version");
    return 0;
}

static int put(const char *name, hid_t ftype, hid_t mtype, const void *data)
{
    int c, bad = 0;
    const column *col = NULL;
    hsize_t dims[1];
    hid_t space, dset;

    for (c = 0; c < NCOLS; c++)
        if (strcmp(COLS[c].name, name) == 0) col = &COLS[c];
    if (!col) return fail("unknown column");
    if (S.particles < 0) return fail("no file is open");

    dims[0] = S.ntot;
    space = H5Screate_simple(1, dims, NULL);
    if (space < 0) return fail("H5Screate_simple");
    dset = H5Dcreate2(S.particles, name, ftype, space, H5P_DEFAULT,
                      H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    if (dset < 0) return fail("H5Dcreate2");
    if (H5Dwrite(dset, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0)
        bad = 1;
    if (!bad && (attr_string(dset, "description", col->description) ||
                 attr_string(dset, "units", col->units)))
        bad = 1;
    H5Dclose(dset);
    return bad ? fail("writing a dataset") : 0;
}

/* The caller passes the first element of a Fortran array of at least ntot. */
int ss_h5_put_double(const char *name, const double *a)
{
    return put(name, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, a);
}

int ss_h5_put_int(const char *name, const int *a)
{
    return put(name, H5T_STD_I32LE, H5T_NATIVE_INT, a);
}

int ss_h5_close(void)
{
    herr_t r = 0;
    if (S.particles >= 0) r |= H5Gclose(S.particles);
    if (S.header >= 0) r |= H5Gclose(S.header);
    if (S.file >= 0) r |= H5Fclose(S.file);
    S.particles = S.header = S.file = -1;
    if (r < 0) {
        remove(S.tmp);
        fprintf(stderr, " hdf5_dump: closing %s failed\n", S.tmp);
        return 1;
    }
    if (rename(S.tmp, S.path) != 0) {
        fprintf(stderr, " hdf5_dump: cannot rename %s to %s\n", S.tmp,
                S.path);
        return 1;
    }
    return 0;
}
