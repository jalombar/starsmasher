#!/usr/bin/env python3
"""StarSmasher test suite.

Runs short calculations and checks them against expectations that must hold,
rather than against stored output.  See tests/README.md.
"""
import argparse, glob, math, os, random, re, shutil, signal, struct, subprocess, sys, tempfile, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(ROOT, 'parallel_bleeding_edge', 'src')

class Fail(Exception):
    pass

class Skip(Exception):
    """The test cannot say anything here (a missing build option or library)."""

# ---------------------------------------------------------------- helpers ---

def sh(cmd, cwd=None, timeout=1800, check=True, kill_group=False, env=None):
    if kill_group:
        # a run that hangs (the ranks wait for each other) is killed with its ranks, not only mpirun, which would leave them spinning
        q = subprocess.Popen(cmd, cwd=cwd, shell=isinstance(cmd, str), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             text=True, start_new_session=True, env=env)
        try:
            out, err = q.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(q.pid, signal.SIGKILL)
            q.communicate()
            raise
        p = subprocess.CompletedProcess(cmd, q.returncode, out, err)
    else:
        p = subprocess.run(cmd, cwd=cwd, shell=isinstance(cmd, str),
                           capture_output=True, text=True, timeout=timeout, env=env)
    if check and p.returncode != 0:
        raise Fail('command failed: %s\n%s' % (cmd, (p.stderr or p.stdout)[-1500:]))
    return p

def build(target='gpu', base=None, hdf5=False):
    """Build in a copy of the source tree, never in the repository itself.

    Building in place would leave object files, .mod files and executables
    scattered through parallel_bleeding_edge/, which is both untidy and a
    trap for anyone who then runs `git add -A`.
    """
    work = os.path.join(base, 'build-%s' % target, 'src')
    shutil.copytree(SRC, work)
    p = sh('make %s%s' % ('' if target == 'gpu' else target,
                          ' HDF5=1' if hdf5 else ''), cwd=work, check=False)
    parent = os.path.dirname(work)
    exe = [os.path.join(parent, f) for f in os.listdir(parent)
           if f.endswith('_%s_sph' % target)]
    if not exe:
        raise Fail('build produced no _%s_sph executable\n%s'
                   % (target, (p.stdout + p.stderr)[-1500:]))
    return exe[0]

def sph_input(cpu=False, **kw):
    # STARSMASHER_TEST_EXTRA="ngravtree=1,treeacc=0.45" adds namelist entries to every run (the
    # settings of a test win), to run the suite with another gravity solver; the reference values
    # of the tests are not changed
    for item in filter(None, os.environ.get('STARSMASHER_TEST_EXTRA', '').split(',')):
        k, v = item.split('=', 1)
        kw.setdefault(k.strip(), v.strip())
    # the CPU build cannot do ngravtree=2 (the visit of the tree on the GPUs): it gets ngravtree=1,
    # which has the same physics
    if cpu and str(kw.get('ngravtree', '')).strip() == '2':
        kw['ngravtree'] = '1'
    body = '\n'.join(' %s=%s,' % (k, v) for k, v in kw.items())
    return ' &input\n%s\n &end\n' % body

def run(exe, workdir, np=4, init="'1es'", timeout=1800, threads=None, **settings):
    # threads: OpenMP threads of each rank for this run (the default is OMP_NUM_THREADS of the environment)
    env = None
    if threads:
        env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    os.makedirs(workdir, exist_ok=True)
    shutil.copy(exe, os.path.join(workdir, 'run_sph'))
    open(os.path.join(workdir, 'sph.init'), 'w').write(
        ' &INITT\n INAME=%s\n &END\n' % init)
    open(os.path.join(workdir, 'sph.input'), 'w').write(
        sph_input(cpu=exe.endswith('_cpu_sph'), **settings))
    p = sh(['mpirun', '--oversubscribe'] + (['--bind-to', 'none'] if threads else []) + ['-np', str(np), './run_sph'],
           cwd=workdir, timeout=timeout, check=False, kill_group=True, env=env)
    return p, os.path.join(workdir, 'energy0.sph')

def energies(path):
    """Rows of energy0.sph as (t, W, T, U, Etot)."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        raise Fail('no energy0.sph was written')
    rows = []
    for line in open(path):
        f = line.split()
        if len(f) >= 5:
            try:
                rows.append(tuple(float(x.replace('D', 'E')) for x in f[:5]))
            except ValueError:
                pass
    if not rows:
        raise Fail('energy0.sph contained no usable rows')
    return rows

def close(a, b, tol, what):
    if b == 0:
        raise Fail('%s: reference value is zero' % what)
    rel = abs(a - b) / abs(b)
    if rel > tol:
        raise Fail('%s: got %.8g, expected %.8g (relative %.2e > %.2e)'
                   % (what, a, b, rel, tol))
    return rel

def energy_words(path):
    """The first five words of each row of energy0.sph as written (t, W, T, U, Etot), for the comparisons that are made in units of the last printed digit."""
    rows = []
    if os.path.exists(path):
        for line in open(path):
            f = line.split()
            if len(f) >= 5:
                rows.append(f[:5])
    return rows

def last_digit(word):
    """The value of the last digit of a number as the Fortran of energy0.sph prints it (7 digits: -0.1714286, 0.8571429E-01): 1e-7 and 1e-8 there."""
    mant, _, ex = word.upper().replace('D', 'E').partition('E')
    dec = len(mant.split('.')[1]) if '.' in mant else 0
    return 10.0 ** (int(ex or 0) - dec)

def close_digits(a, b, units, what):
    """a and b are numbers as written in energy0.sph: they must differ by at most `units` units of the last printed digit (of the larger one)."""
    d = abs(float(a.upper().replace('D', 'E')) - float(b.upper().replace('D', 'E'))) / max(last_digit(a), last_digit(b))
    if d > units * (1.0 + 1e-9):
        raise Fail('%s: %s and %s differ by %.2f units of the last printed digit (limit %g)' % (what, a, b, d, units))
    return d

# ------------------------------------------------------------------ tests ---
# A polytrope small enough to finish quickly but large enough to be meaningful.
POLY = dict(n=10000, tf=10, dtout=1, nrelax=1, treloff=0,
            starmass='0.2d0', starradius='0.2d0', ngravprocs=1)

def t_virial(ctx):
    """A relaxed n=1.5 polytrope has W = -3/(5-n) GM^2/R and U = -W/2."""
    _, e = run(ctx.exe, ctx.dir('virial'), np=ctx.np, **POLY)
    t, W, T, U, Etot = energies(e)[0]
    M = R = 0.2
    expected_W = -3.0 / (5.0 - 1.5) * M * M / R
    close(W, expected_W, 1e-3, 'gravitational potential energy')
    close(U, -W / 2.0, 1e-3, 'internal energy (virial theorem)')
    return 'W=%.6f (theory %.6f), U=%.6f (-W/2=%.6f)' % (W, expected_W, U, -W/2)

def t_energy(ctx):
    """Total energy is conserved through the run."""
    _, e = run(ctx.exe, ctx.dir('energy'), np=ctx.np, **POLY)
    rows = energies(e)
    first, last = rows[0][4], rows[-1][4]
    rel = close(last, first, 1e-3, 'total energy drift')
    return 'Etot %.8g -> %.8g over %d rows (drift %.2e)' % (first, last, len(rows), rel)

def t_quiet(ctx):
    """A relaxed star stays put: kinetic energy stays negligible."""
    _, e = run(ctx.exe, ctx.dir('quiet'), np=ctx.np, **POLY)
    rows = energies(e)
    T, U = rows[-1][2], rows[-1][3]
    if abs(T) > 1e-3 * abs(U):
        raise Fail('kinetic energy %.3e is not negligible against internal %.3e' % (T, U))
    return 'final T=%.3e, which is %.1e of U' % (T, abs(T / U))

def t_ranks(ctx):
    """The answer must not depend on the number of MPI ranks."""
    out = {}
    for np_ in (1, ctx.np):
        _, e = run(ctx.exe, ctx.dir('ranks%d' % np_), np=np_, **POLY)
        out[np_] = energies(e)[0]
    a, b = out[1], out[ctx.np]
    for i, name in ((1, 'W'), (3, 'U'), (4, 'Etot')):
        close(b[i], a[i], 1e-6, '%s at np=%d vs np=1' % (name, ctx.np))
    return 'np=1 and np=%d agree on W, U and Etot to 1e-6' % ctx.np

def tree2_extra():
    """True when STARSMASHER_TEST_EXTRA asks for ngravtree=2 or -1: the GPU build then runs the tree on the GPU
    (with -1, below ngravtreeauto, the direct sum) and the CPU build the tree on the CPUs (ngravtree=1), which
    differ by the rounding of the float arithmetic (and by the error of the tree).  The stars that the two builds
    make in the relaxation differ by that, amplified by the opening of the nodes (a tree force is not continuous
    in the positions), so cpu_gpu and block_cpu_gpu compare them with wider limits: 5e-6 and 2e-4, against 1e-6
    and 1e-5 for the direct sum.  The worst differences measured in sixteen runs (the star of 2000 and of 10000
    particles, 2 and 4 ranks, with and without OpenMP, ngravtree=2 and -1) are 1.2e-6 in cpu_gpu and 5.6e-5 in
    block_cpu_gpu."""
    e = os.environ.get('STARSMASHER_TEST_EXTRA', '').replace(' ', '')
    return 'ngravtree=2' in e or 'ngravtree=-1' in e


def t_cpu_gpu(ctx):
    """The CPU and GPU builds must agree."""
    if not ctx.cpu_exe:
        raise Fail('no CPU build available (pass --with-cpu to build one)')
    res = {}
    for tag, exe in (('gpu', ctx.exe), ('cpu', ctx.cpu_exe)):
        _, e = run(exe, ctx.dir('cpugpu_' + tag), np=ctx.np, **POLY)
        res[tag] = energies(e)[0]
    tol = 5e-6 if tree2_extra() else 1e-6
    for i, name in ((1, 'W'), (3, 'U'), (4, 'Etot')):
        close(res['cpu'][i], res['gpu'][i], tol, '%s, cpu vs gpu' % name)
    return 'cpu and gpu agree on W, U and Etot to %.0e' % tol

def t_block_cpu_gpu(ctx):
    """With block timesteps the CPU and GPU builds must agree.

    An unrelaxed polytrope run dynamically puts its particles on several
    timestep bins, so most substeps have only some particles active.  The GPU
    build then computes gravity on the active particles only, and the CPU
    build does the same with its own active-set sum.
    """
    if not ctx.cpu_exe:
        raise Fail('no CPU build available (pass --with-cpu to build one)')
    s = dict(POLY, nrelax=0, nblock=1, tf=1, dtout=0.5)
    res = {}
    for tag, exe in (('gpu', ctx.exe), ('cpu', ctx.cpu_exe)):
        d = ctx.dir('block_cpugpu_' + tag)
        p, e = run(exe, d, np=ctx.np, **s)
        log = os.path.join(d, 'log0.sph')
        text = open(log).read() if os.path.exists(log) else ''
        active = [float(x) for x in re.findall(r'mean active=\s*(\S+)', text)]
        hist = re.findall(r'bin histogram \(bins 0\.\.finest\):([\d\s]+)', text)
        if not active or not hist:
            raise Fail('the %s run wrote no block sync lines\n' % tag + p.stdout[-800:])
        ntot = sum(int(x) for x in hist[0].split())
        if min(active) >= ntot:
            raise Fail('every particle was active in every substep (%s build), so '
                       'the active-particle gravity was never used' % tag)
        res[tag] = energies(e), min(active) / ntot
    rows, frac = res['gpu']
    close(rows[-1][4], rows[0][4], 1e-3, 'total energy with block steps (gpu)')
    tol = 2e-4 if tree2_extra() else 1e-5
    for i, name in ((1, 'W'), (3, 'U'), (4, 'Etot')):
        close(res['cpu'][0][-1][i], rows[-1][i], tol, '%s at t=%g, cpu vs gpu' % (name, rows[-1][0]))
    return ('t=%g: cpu and gpu agree on W, U and Etot to %.0e '
            '(as few as %.0f%% of particles active per substep)' % (rows[-1][0], tol, 100 * frac))

def t_gpu_hilo(ctx):
    """GPU gravity with ngravhilo=1 (the default) and 0 agrees for one star.

    A single star sits on the reference point, where the float remainders of
    the positions change the forces only in the last digits, so both settings
    must give the same answer.  This keeps ngravhilo=0 from going stale.
    """
    res = {}
    for flag in (1, 0):
        _, e = run(ctx.exe, ctx.dir('hilo%d' % flag), np=ctx.np, ngravhilo=flag, **POLY)
        res[flag] = energies(e)[0]
    for i, name in ((1, 'W'), (3, 'U'), (4, 'Etot')):
        close(res[0][i], res[1][i], 1e-6, '%s, ngravhilo=0 vs 1' % name)
    return 'ngravhilo=0 and 1 agree on W, U and Etot to 1e-6'

def t_gpu_tree(ctx):
    """The tree gravity on the GPU (ngravtree=2) agrees with the tree on the CPUs (ngravtree=1).

    For the same opening angle and leaf size the two have the same physics (the criterion, the
    Taylor series of the far field, the pair formula), so the accelerations, the potentials and the
    energies of the same stars differ by rounding only: the GPU does the arithmetic in float.  They
    are compared on the same state: a star is relaxed once, and each solver starts from that
    snapshot (as the test startfiles does with tri).  Two stars made in two runs would not do: a
    tree force is not continuous in the positions (a node is opened or not), so two runs that
    differ by rounding follow different paths, and the stars they make differ by the error of the
    tree.
    """
    d = ctx.dir('gputree_relax')
    s = dict(POLY); s['tf'] = 2; s['dtout'] = 1; s['n'] = 2000
    p, e = run(ctx.exe, d, np=ctx.np, **s)
    snaps = sorted(glob.glob(os.path.join(d, 'out*.sph')))
    if not snaps:
        raise Fail('the relaxation produced no snapshot to use as a start file\n'
                   + (p.stdout + p.stderr)[-800:])
    m = POLY['starmass'].replace('d0', '')
    res = {}
    # gpuh: the same on the GPU with a threshold so low that most leaves take the path of the leaves
    # with a very long list (divided in chunks), see grav_tree_gpu.cu
    for tag, nt in (('cpu', 1), ('gpu', 2), ('gpuh', 2)):
        d2 = ctx.dir('gputree_' + tag)
        for name in ('body_one.sph', 'body_two.sph', 'body_three.sph'):
            shutil.copy(snaps[-1], os.path.join(d2, name))
        open(os.path.join(d2, 'three_bodies.txt'), 'w').write(
            'triple:\n%s\n0.0 0.0 0.0\n0.0 0.0 0.0\n%s\n20.0 0.0 0.0\n0.0 0.3 0.0\n'
            '%s\n0.0 60.0 0.0\n-0.2 0.0 0.0\n' % (m, m, m))
        if tag == 'gpuh':
            os.environ['STARSMASHER_TREE_HEAVY'] = '24'
        try:
            p2, e2 = run(ctx.exe, d2, np=ctx.np, init="'tri'", tf=0.05, dtout=0.05, nrelax=0,
                         startfile1="'body_one.sph'", startfile2="'body_two.sph'",
                         startfile3="'body_three.sph'", triplefile="'three_bodies.txt'",
                         ngravtree=nt, treeacc='0.4d0', nleafgrav=16, ngravprocs=1)
        finally:
            os.environ.pop('STARSMASHER_TREE_HEAVY', None)
        dump = os.path.join(d2, 'out0000.sph')
        if not os.path.exists(dump):
            raise Fail('the ngravtree=%d run wrote no dump\n' % nt + (p2.stdout + p2.stderr)[-800:])
        res[tag] = energies(e2), read_dump(dump, full=True)[1]
    ec, dc = res['cpu']
    out = []
    for tag in ('gpu', 'gpuh'):
        eg, dg = res[tag]
        for i, name in ((1, 'W'), (3, 'U'), (4, 'Etot')):
            close(eg[0][i], ec[0][i], 1e-6, '%s at t=0, ngravtree=2 vs 1 (%s)' % (name, tag))
        # the two start from the same snapshot, so the positions are the same at the first dump
        # to the change of one step
        dx = max(abs(a[c] - b[c]) for a, b in zip(dc, dg) for c in (0, 1, 2))
        if dx > 1e-6:
            raise Fail('the two runs do not start from the same state (positions differ by %.2e)'
                       % dx)
        # the accelerations are compared with the scale of the gravity of a star, G M / R^2
        g0 = (float(POLY['starmass'].replace('d0', ''))
              / float(POLY['starradius'].replace('d0', '')) ** 2)
        n = len(dc)
        sq = mx = 0.0
        spot = 0.0
        for a, b in zip(dc, dg):
            d2s = (a[9] - b[9]) ** 2 + (a[10] - b[10]) ** 2 + (a[11] - b[11]) ** 2
            sq += d2s
            mx = max(mx, math.sqrt(d2s))
            spot += ((a[14] - b[14]) / a[14]) ** 2
        rms = math.sqrt(sq / n) / g0
        prms = math.sqrt(spot / n)
        # the limits are three times the worst difference measured in eight runs (2 and 4 ranks,
        # with and without OpenMP, ngravtree=2 and -1 in the environment): accelerations rms 9e-7,
        # max 5.5e-6 of G M / R^2, potentials 1.9e-8, energies identical
        if rms > 5e-6 or mx / g0 > 2e-5:
            raise Fail('accelerations differ (%s): rms %.2e, max %.2e of G M / R^2 '
                       '(limits 5e-6, 2e-5)' % (tag, rms, mx / g0))
        if prms > 1e-7:
            raise Fail('potentials differ (%s): rms relative %.2e (limit 1e-7)' % (tag, prms))
        out.append('%s: accelerations rms %.1e max %.1e of G M/R^2, potentials rms %.1e'
                   % (tag, rms, mx / g0, prms))
    return ('ngravtree=2 and 1 agree on the same state of %d particles (W, U, Etot to 1e-6), '
            'also with the leaves divided in chunks (gpuh) - ' % n) + '; '.join(out)

def t_guard_ngravprocs(ctx):
    """Requesting more gravity processes than GPUs must be clamped, not fatal."""
    s = dict(POLY); s['ngravprocs'] = -8; s['tf'] = 2
    p, e = run(ctx.exe, ctx.dir('ngrav'), np=ctx.np, **s)
    if not os.path.exists(e) or os.path.getsize(e) == 0:
        raise Fail('run with ngravprocs=-8 produced no output; it should have been '
                   'clamped to the number of GPUs present\n' + p.stdout[-800:])
    energies(e)
    return 'ngravprocs=-8 was clamped and the run completed'


def t_startfiles(ctx):
    """Bodies must be read from the names sph.input gives, not from fixed ones.

    Several setup routines used to open hardcoded file names, which meant a
    directory could only ever hold one calculation.  This relaxes a star, hands
    the same snapshot to `tri` three times under names that are deliberately not
    the defaults, and checks all three were read.  If a routine falls back to a
    hardcoded name the run cannot start, because no such file exists here.
    """
    import glob
    d = ctx.dir('startfiles')
    s = dict(POLY); s['tf'] = 2; s['dtout'] = 1; s['n'] = 2000
    p, e = run(ctx.exe, d, np=ctx.np, **s)
    snaps = sorted(glob.glob(os.path.join(d, 'out*.sph')))
    if not snaps:
        raise Fail('the relaxation produced no snapshot to use as a start file\n'
                   + p.stdout[-800:])
    # A fresh directory for the second run.  The relaxation leaves a
    # restartrad.sph behind, and the code picks that up automatically, so
    # reusing the directory would silently resume the relaxation instead of
    # running tri at all.
    d2 = ctx.dir('startfiles_tri')
    os.makedirs(d2, exist_ok=True)
    for name in ('body_one.sph', 'body_two.sph', 'body_three.sph'):
        shutil.copy(snaps[-1], os.path.join(d2, name))

    # tri checks each mass against the start file to one part in 1e8, and all
    # three bodies here are the same relaxed POLY star.
    m = POLY['starmass'].replace('d0', '')
    open(os.path.join(d2, 'three_bodies.txt'), 'w').write(
        'triple:\n%s\n0.0 0.0 0.0\n0.0 0.0 0.0\n'
        '%s\n20.0 0.0 0.0\n0.0 0.3 0.0\n'
        '%s\n0.0 60.0 0.0\n-0.2 0.0 0.0\n' % (m, m, m))

    p2, _ = run(ctx.exe, d2, np=ctx.np, init="'tri'",
                tf=0.1, dtout=0.1, nrelax=0,
                startfile1="'body_one.sph'",
                startfile2="'body_two.sph'",
                startfile3="'body_three.sph'",
                triplefile="'three_bodies.txt'")
    log = os.path.join(d2, 'log0.sph')
    text = open(log).read() if os.path.exists(log) else p2.stdout
    m = re.search(r'triple:\s*n=\s*(\d+)\s*ntot=\s*(\d+)', text)
    if not m:
        raise Fail('tri did not report a particle count; it probably could not '
                   'open one of the start files\n' + text[-800:])
    ntot = int(m.group(2))
    if ntot % 3 != 0:
        raise Fail('tri read %d particles, which is not three equal bodies' % ntot)
    return 'tri read three bodies (%d particles) from non-default file names' % ntot

# ------------------------------------------------------- periodic boxes ---
#
# A periodic run starts from particles written by a script (iname='pbx'), in
# the out*.sph layout: a header record, one record per particle, and ntot
# again.  These helpers write and read that layout with nothing but struct.

HEADER_FMT = '<2i5d2id i3d2i3di d i3d'
PART_FMT = '<16d i d'

def write_dump(path, parts):
    """parts: list of dicts with x,y,z,am,hp,u and optionally vx,vy,vz."""
    def rec(f, b):
        f.write(struct.pack('<i', len(b)) + b + struct.pack('<i', len(b)))
    hdr = struct.pack(HEADER_FMT, len(parts), 23, 0, 0, 0, 0, 0, 0, 0, 0,
                      3, 1, 2, 1e30, 0, 0, 1e30, 0, 0, 0, 0, 0, 0, 0, 0)
    with open(path, 'wb') as f:
        rec(f, hdr)
        for p in parts:
            rec(f, struct.pack(PART_FMT, p['x'], p['y'], p['z'], p['am'], p['hp'],
                               1.0, p.get('vx', 0), p.get('vy', 0), p.get('vz', 0),
                               0, 0, 0, p['u'], 0, 0, 1.0, 0, 0))
        rec(f, struct.pack('<i', len(parts)))

def read_dump(path, full=False):
    """(t, list of (x, y, z, am, rho, vx, vy, vz, u)); with full, the first 16 doubles of every
    particle: x y z am hp rho vx vy vz vxdot vydot vzdot u udot grpot meanmolecular"""
    data = open(path, 'rb').read()
    n = struct.unpack_from('<i', data, 4)[0]
    t = struct.unpack_from('<d', data, 4 + 4 * 4 + 5 * 8)[0]
    pos = 4 + struct.unpack_from('<i', data, 0)[0] + 4
    out = []
    for _ in range(n):
        v = struct.unpack_from(PART_FMT, data, pos + 4)
        out.append(v[:16] if full else (v[0], v[1], v[2], v[3], v[5], v[6], v[7], v[8], v[12]))
        pos += struct.calcsize(PART_FMT) + 8
    return t, out

def hcp_box(nx, ny, nz, jitter, seed):
    """hcp lattice of unit spacing in a box commensurate with it (ny, nz
    even), centred on the origin, slightly jittered."""
    rng = random.Random(seed)
    L = (nx, ny * math.sqrt(3) / 2, nz * math.sqrt(2.0 / 3.0))
    pts = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                p = [i + 0.5 * (j % 2) + 0.5 * (k % 2),
                     (j + (k % 2) / 3.0) * math.sqrt(3) / 2,
                     k * math.sqrt(2.0 / 3.0)]
                pts.append([(p[d] + jitter * rng.uniform(-1, 1)) % L[d] - L[d] / 2
                            for d in range(3)])
    return pts, L

PBOX = dict(ngr=0, neos=0, gam='1.6666666666666667d0', nintvar=2, nrelax=0,
            startfile1="'box.sph'")

def t_periodic(ctx):
    """A periodic box conserves momentum and energy and has no seams at its faces."""
    d = ctx.dir('periodic')
    pts, L = hcp_box(10, 10, 10, 0.1, 1)
    rng = random.Random(2)
    n = len(pts)
    vel = [[rng.gauss(0, 0.05) for _ in range(3)] for _ in pts]
    mean = [sum(v[c] for v in vel) / n for c in range(3)]
    parts = [dict(x=p[0], y=p[1], z=p[2], am=L[0] * L[1] * L[2] / n, hp=1.0, u=1.0,
                  vx=v[0] - mean[0], vy=v[1] - mean[1], vz=v[2] - mean[2])
             for p, v in zip(pts, vel)]
    write_dump(os.path.join(d, 'box.sph'), parts)
    # cn1 is small so that the leapfrog's own truncation error, which with
    # the default takes ~5 steps across this box, does not mask the check
    s = dict(PBOX, tf=1, dtout=1, cn1=0.05, boxx=L[0], boxy=L[1], boxz=L[2])
    p, e = run(ctx.exe, d, np=ctx.np, init="'pbx'", **s)
    dumps = sorted(f for f in os.listdir(d) if f.startswith('out') and f.endswith('.sph'))
    if len(dumps) < 2:
        raise Fail('the periodic run wrote too few dumps\n' + p.stdout[-800:])
    t0, first = read_dump(os.path.join(d, dumps[0]))
    t1, last = read_dump(os.path.join(d, dumps[-1]))
    # momentum: every pair force is antisymmetric, so the total stays zero
    # to round-off however the particles move across the faces
    mtot = sum(q[3] for q in last)
    vrms = math.sqrt(sum(q[3] * (q[5] ** 2 + q[6] ** 2 + q[7] ** 2) for q in last) / mtot)
    pmax = max(abs(sum(q[3] * q[5 + c] for q in last)) for c in range(3))
    if pmax > 1e-10 * mtot * vrms:
        raise Fail('total momentum drifted to %.3e (scale %.3e)' % (pmax, mtot * vrms))
    rows = energies(e)
    dE = abs(rows[-1][4] - rows[0][4]) / abs(rows[0][4])
    close(rows[-1][4], rows[0][4], 1e-3, 'total energy')
    for q in last:
        if any(abs(q[c]) > L[c] / 2 * (1 + 1e-12) for c in range(3)):
            raise Fail('a particle left the box: %r' % (q[:3],))
    # a seam would show as a density step for particles near a face
    rho_face = [q[4] for q in first if any(abs(q[c]) > L[c] / 2 - 1 for c in range(3))]
    rho_in = [q[4] for q in first if all(abs(q[c]) < L[c] / 2 - 1 for c in range(3))]
    close(sum(rho_face) / len(rho_face), sum(rho_in) / len(rho_in), 0.01,
          'density near the faces against the interior')
    return ('t=%.2f: momentum %.1e of scale, dE/E=%.1e, face/interior density %.4f'
            % (t1, pmax / (mtot * vrms), dE,
               sum(rho_face) / len(rho_face) / (sum(rho_in) / len(rho_in))))

def t_periodic_block(ctx):
    """Block timesteps in a periodic box agree with shared ones and keep every particle in the box."""
    pts, L = hcp_box(10, 10, 10, 0.1, 4)
    n = len(pts)
    # a hot half drives a flow across the faces and puts the two halves on
    # different timestep bins, so that block steps really are mixed
    parts = [dict(x=p[0], y=p[1], z=p[2], am=L[0] * L[1] * L[2] / n, hp=1.0,
                  u=4.0 if p[0] > 0 else 1.0) for p in pts]
    out = {}
    for tag, nblock in (('shared', 0), ('block', 1)):
        d = ctx.dir('periodic_block_' + tag)
        write_dump(os.path.join(d, 'box.sph'), parts)
        s = dict(PBOX, tf=0.5, dtout=0.5, nblock=nblock, boxx=L[0], boxy=L[1], boxz=L[2])
        p, e = run(ctx.exe, d, np=ctx.np, init="'pbx'", **s)
        dumps = sorted(f for f in os.listdir(d) if f.startswith('out') and f.endswith('.sph'))
        if len(dumps) < 2:
            raise Fail('the %s run wrote too few dumps\n' % tag + p.stdout[-800:])
        out[tag] = read_dump(os.path.join(d, dumps[-1])), energies(e)
    (tb, last), rows = out['block']
    for q in last:
        if any(abs(q[c]) > L[c] / 2 * (1 + 1e-12) for c in range(3)):
            raise Fail('a particle left the box: %r' % (q[:3],))
    close(rows[-1][4], rows[0][4], 1e-2, 'total energy with block steps')
    (ts, ref), _ = out['shared']
    drho = sum(abs(a[4] - b[4]) / b[4] for a, b in zip(last, ref)) / n
    if drho > 0.01:
        raise Fail('block and shared steps differ in density by %.2e on average (t=%.3f, %.3f)'
                   % (drho, tb, ts))
    return ('t=%.2f: mean |rho_block-rho_shared|/rho = %.1e, dE/E=%.1e'
            % (tb, drho, abs(rows[-1][4] - rows[0][4]) / abs(rows[0][4])))

def t_guard_periodic(ctx):
    """A periodic box too small for the kernel stops with a message."""
    d = ctx.dir('guard_periodic')
    pts, L = hcp_box(4, 4, 4, 0.0, 3)
    parts = [dict(x=p[0], y=p[1], z=p[2], am=L[0] * L[1] * L[2] / len(pts), hp=1.0, u=1.0)
             for p in pts]
    write_dump(os.path.join(d, 'box.sph'), parts)
    s = dict(PBOX, tf=0.1, dtout=0.1, boxx=L[0], boxy=L[1], boxz=L[2])
    p, _ = run(ctx.exe, d, np=ctx.np, init="'pbx'", **s)
    text = p.stdout + p.stderr
    if 'exceeds periodic padding' not in text:
        raise Fail('expected the padding guard to stop the run\n' + text[-800:])
    return 'stopped with: search radius exceeds periodic padding'

# The automatic values (ngravtree=-1, ngpusph=-1).  The rule is in init.f (auto_limits, gravtree_auto, sph_auto) and in docs/source/using/gravity.rst; the limits, from results/g_regola.csv of the phase 10
# (relaxed Suns of 9966 to 499860 particles, 7 to 63 threads per rank, one rank with one GPU): the smallest star measured is 9966 particles, below it the choice of the phases 7 and 8 stays (the direct
# sum; ngpusph 7 with at most 8 threads per rank, 0 with 32 or more, in between 7 below 300000).  From there up, with one rank and one GPU:
#   ngravtree=-1 (ngravtreeauto=-1): the tree on the GPU from 9966 particles with ngpusph=-1 (the pair of the tree and ngpusph=7 is the fastest) or 3 or 7 written, from 25000 with ngpusph=1, from 75000 with ngpusph=0 and
#   48 threads or more (450000 with fewer); two GPUs (2 ranks): from 15000 with the SPH on the GPU, from 450000 with ngpusph 0 or 1; more GPUs (4 ranks here): from 299000
#   ngpusph=-1: with the tree 7; the direct sum on one GPU 7 with at most 8 threads, with 9 to 31 threads below 450000 particles, with 32 to 47 below 150000 and with 48 or more below 75000 (0 from there); on two GPUs
#   7 (1 from 450000); on more 1
# (the rank counts of the cases do not depend on the --np of the suite.)  Each case is run with the automatic values, the pair that was chosen is read from log0.sph (and compared with the rule above), and the run is made again with that pair written by hand: the energies of
# every row and the bytes of every dump must be those of the automatic run.  The stars are polytropes of n particles (a few less come out of the lattice), a factor of 1.1 or more from the limit.
AUTO_CASES_GRAVITY = [   # (label, ranks, threads, n, written, wanted (ngravtree, ngpusph))
    ('both_small_1t',   1, 1, 5000,  {}, (0, 7)),
    ('both_big_1t',     1, 1, 14000, {}, (2, 7)),
    ('both_small_63t',  1, 63, 5000, {}, (0, 0)),
    ('both_big_63t',    1, 63, 14000, {}, (2, 7)),
    ('s0_75k_below',    1, 63, 60000, dict(ngpusph=0), (0, 0)),
    ('s0_75k_above',    1, 63, 90000, dict(ngpusph=0), (2, 0)),
    ('s0_1t_big',       1, 1, 14000, dict(ngpusph=0), (0, 0)),
    ('s1_below',        1, 1, 18000, dict(ngpusph=1), (0, 1)),
    ('s1_above',        1, 1, 33000, dict(ngpusph=1), (2, 1)),
    ('s7_small',        1, 1, 5000,  dict(ngpusph=7), (0, 7)),
    ('s7_big',          1, 1, 14000, dict(ngpusph=7), (2, 7)),
    ('s3_big',          1, 1, 14000, dict(ngpusph=3), (2, 3)),
    ('two_ranks_small', 2, 1, 12000, {}, (0, 7)),
    ('two_ranks_big',   2, 1, 18000, {}, (2, 7)),
    ('two_ranks_s7',    2, 1, 18000, dict(ngpusph=7), (2, 7)),
    ('two_ranks_s1',    2, 1, 33000, dict(ngpusph=1), (0, 1)),
    ('two_ranks_s0',    2, 1, 33000, dict(ngpusph=0), (0, 0)),
    ('four_ranks',      4, 1, 14000, {}, (0, 1)),
    ('nleafgrav_40',    1, 1, 14000, dict(nleafgrav=40), (0, 7)),
    ('no_gravity',      1, 1, 14000, dict(ngr=0), (0, 7)),
    ('user_below',      1, 1, 14000, dict(ngravtreeauto=20000), (0, 7)),
    ('user_above',      1, 1, 14000, dict(ngravtreeauto=10000), (2, 7)),
]
AUTO_CASES_SPH = [       # (label, ranks, threads, n, written, wanted ngpusph)
    ('direct_1t',         1, 1, 14000, dict(ngravtree=0), 7),
    ('direct_16t',        1, 16, 14000, dict(ngravtree=0), 7),
    ('direct_32t_small',  1, 32, 5000, dict(ngravtree=0), 0),
    ('direct_32t',        1, 32, 14000, dict(ngravtree=0), 7),
    ('direct_32t_below',  1, 32, 136000, dict(ngravtree=0), 7),
    ('direct_32t_above',  1, 32, 166000, dict(ngravtree=0), 0),
    ('direct_63t_below',  1, 63, 68000, dict(ngravtree=0), 7),
    ('direct_63t_above',  1, 63, 83000, dict(ngravtree=0), 0),
    ('direct_63t_small',  1, 63, 5000, dict(ngravtree=0), 0),
    ('direct_2_ranks',    2, 1, 14000, dict(ngravtree=0), 7),
    ('direct_2_ranks_small', 2, 1, 5000, dict(ngravtree=0), 1),
    ('direct_4_ranks',    4, 1, 14000, dict(ngravtree=0), 1),
    ('tree_gpu',          1, 1, 14000, dict(ngravtree=2, treeacc='0.5d0'), 7),
    ('tree_cpu',          1, 1, 14000, dict(ngravtree=1, treeacc='0.5d0'), 7),
    ('tree_4_ranks',      4, 1, 14000, dict(ngravtree=1, treeacc='0.5d0'), 7),
]

def periodic_box_setup(d):
    """The dump of a small periodic box (as t_periodic) in directory d; returns the settings of the run."""
    pts, L = hcp_box(10, 10, 10, 0.1, 1)
    n = len(pts)
    parts = [dict(x=q[0], y=q[1], z=q[2], am=L[0] * L[1] * L[2] / n, hp=1.0, u=1.0, vx=0.0, vy=0.0, vz=0.0) for q in pts]
    write_dump(os.path.join(d, 'box.sph'), parts)
    return dict(PBOX, tf=0.2, dtout=1, cn1=0.05, boxx=L[0], boxy=L[1], boxz=L[2])

def auto_read(d):
    """What the automatic values chose, read from log0.sph: (ngravtree, ngpusph, the text); a value that was written (not -1) is None."""
    log = os.path.join(d, 'log0.sph')
    text = open(log).read() if os.path.exists(log) else ''
    mg = re.search(r'gravity: ngravtree=-1 \(auto\):[^\n]*-> ngravtree=(\d)', text)
    ms = re.search(r'SPH on the GPU: ngpusph=-1 \(auto\):[^\n]*-> ngpusph=(\d+)', text)
    return (int(mg.group(1)) if mg else None), (int(ms.group(1)) if ms else None), text

def auto_case(ctx, exe, tag, np_, threads, n, written, wanted, base=None, init="'1es'", prep=None):
    """One configuration: the run with the automatic values (ngravtree and ngpusph -1 unless `written` gives them), the pair that it chose read from log0.sph and compared with `wanted` (an ngravtree, ngpusph
    pair, or one of them, or None), and then the same run with that pair written by hand: the energies of every row and the bytes of every dump must be the same.  Returns the pair, or None for a case that needs more threads than the build has (no OpenMP)."""
    settings = dict(ngravtree=-1, ngpusph=-1)
    settings.update(written)
    if base is None:
        base = dict(POLY, n=n, nrelax=0, tf=0.2, dtout=1, ngravprocs=1)
    res = {}
    chosen = None
    for mode in ('auto', 'hand'):
        st = dict(base, **settings)
        if mode == 'hand':
            st.update(ngravtree=chosen[0], ngpusph=chosen[1])
        d = ctx.dir('au_%s_%s' % (tag, mode))
        if prep:
            prep(d)
        p, e = run(exe, d, np=np_, threads=threads, init=init, **st)
        g, s, text = auto_read(d)
        if mode == 'auto':
            # a build without OpenMP (make OPENMP=0) has one thread per rank whatever OMP_NUM_THREADS says: a case that needs more threads cannot be tried there
            mt = re.search(r'OpenMP: P=\s*\d+ T=\s*(\d+)', text)
            if threads and threads > 1 and (int(mt.group(1)) if mt else 1) != threads:
                return None
            if (settings['ngravtree'] == -1 and g is None) or (settings['ngpusph'] == -1 and s is None):
                raise Fail('%s: log0.sph does not say what the automatic values chose\n%s' % (tag, (p.stdout + p.stderr)[-600:]))
            chosen = (g if g is not None else int(settings['ngravtree']), s if s is not None else int(settings['ngpusph']))
            if wanted is not None:
                want = tuple(w if w is not None else c for w, c in zip(wanted, chosen))
                if chosen != want:
                    line = [l for l in text.split('\n') if 'auto' in l and ('ngravtree=-1' in l or 'ngpusph=-1' in l)]
                    raise Fail('%s (%d ranks, %d threads, n=%d, written %s): chose ngravtree=%d ngpusph=%d, the rule says %d %d\n%s'
                               % (tag, np_, threads, n, written, chosen[0], chosen[1], want[0], want[1], '\n'.join(line)))
        names = sorted(f for f in os.listdir(d) if f.startswith('out') and f.endswith('.sph'))
        res[mode] = (energies(e), [(f, open(os.path.join(d, f), 'rb').read()) for f in names])
    if res['auto'] != res['hand']:
        raise Fail('%s: ngravtree=-1 / ngpusph=-1 do not give the bits of ngravtree=%d ngpusph=%d written by hand' % (tag, chosen[0], chosen[1]))
    return chosen

def t_gravity_auto(ctx):
    """ngravtree=-1 chooses the gravity from the particles, the threads, the ranks and the SPH, says so in log0.sph, and gives exactly the energies and dumps of the value it chose, set by hand.

    GPU build: with ngravtreeauto written the direct sum below it and the tree walked on the GPU from it up; with the default (ngravtreeauto=-1) the limits measured for the SPH that goes with the gravity
    (both automatic, ngpusph 0, 1, 3 and 7 written), one star below and one above each limit, with one rank and with several; CPU build: the tree on the CPUs.  Each choice is read from log0.sph, compared
    with the rule, and the run is made again with that choice written by hand: the energies of every row and the bytes of every dump must be the same.
    """
    out = []
    n = None
    for tag, exe in (('gpu', ctx.exe), ('cpu', ctx.cpu_exe)):
        if not exe:
            continue
        cpu = exe.endswith('_cpu_sph')
        if n is None:
            # the star has a few particles less than POLY['n'] (the lattice): read the number the code reports
            d0 = ctx.dir('auto_n')
            run(exe, d0, np=ctx.np, ngravtree=-1, ngravtreeauto=10 ** 9, **dict(POLY, tf=0.1, dtout=1))
            m0 = re.search(r'ngravtree=-1 \(auto\):[^\n]*N=(\d+)', open(os.path.join(d0, 'log0.sph')).read())
            if not m0:
                raise Fail('no "gravity: ngravtree=-1 (auto): N=..." line in log0.sph')
            n = int(m0.group(1))
        # ngravtreeauto written: the threshold is the user's, below it the direct sum, at it the tree
        for thr, want_gpu in ((n + 1, 0), (n, 2)):
            want = 1 if cpu else want_gpu
            s = dict(POLY, tf=1, dtout=1)
            d = ctx.dir('auto_%s_%d' % (tag, thr))
            p, e = run(exe, d, np=ctx.np, ngravtree=-1, ngravtreeauto=thr, **s)
            log = os.path.join(d, 'log0.sph')
            text = open(log).read() if os.path.exists(log) else ''
            m = re.search(r'gravity: ngravtree=-1 \(auto\):[^\n]*-> ngravtree=(\d)', text)
            if not m:
                raise Fail('%s build, ngravtreeauto=%d: no "gravity: ngravtree=-1 (auto)" line in '
                           'log0.sph\n' % (tag, thr) + (p.stdout + p.stderr)[-800:])
            if int(m.group(1)) != want:
                raise Fail('%s build, N=%d, ngravtreeauto=%d: chose ngravtree=%s, expected %d\n%s'
                           % (tag, n, thr, m.group(1), want, m.group(0)))
            if 'N=%d' % n not in m.group(0) and not cpu:
                raise Fail('the log line does not say N=%d: %s' % (n, m.group(0)))
            auto = energies(e)
            d2 = ctx.dir('hand_%s_%d' % (tag, thr))
            _, e2 = run(exe, d2, np=ctx.np, ngravtree=want, ngravtreeauto=thr, **s)
            if energies(e2) != auto:
                raise Fail('%s build, ngravtree=-1 and ngravtree=%d by hand do not give the same '
                           'energies' % (tag, want))
            out.append('%s N=%d thr=%d -> %d' % (tag, n, thr, want))
        # the default (ngravtreeauto=-1): the CPU build takes the tree on the CPUs, the GPU build the limits measured
        if cpu:
            chosen = auto_case(ctx, exe, 'cpu_default', 1, 1, 5000, dict(ngpusph=0), (1, None))
            out.append('cpu default -> %d' % chosen[0])
            continue
        for label, np_, threads, nn, written, wanted in AUTO_CASES_GRAVITY:
            chosen = auto_case(ctx, exe, 'g_' + label, np_, threads, nn, written, wanted)
            if chosen is None:
                out.append('%s skipped (no OpenMP)' % label)
                continue
            out.append('%s -> %d,%d' % (label, chosen[0], chosen[1]))
        # a periodic box without gravity: the direct sum, and the SPH on the CPUs (the box is not done on the GPU)
        d0 = ctx.dir('au_pbox_setup')
        chosen = auto_case(ctx, exe, 'g_periodic_box', 1, 1, 0, {}, (0, 0), base=periodic_box_setup(d0), init="'pbx'", prep=periodic_box_setup)
        out.append('periodic_box -> %d,%d' % chosen)
    return 'chosen as the rule says and the same bits as by hand: ' + '; '.join(out)

def t_gpu_sph(ctx):
    """ngpusph=1 (the SPH on the GPU) gives exactly the bits of ngpusph=0, and ngpusphcheck=1 finds no difference.

    With bit 1 of ngpusph the smoothing lengths, the densities and the neighbour lists are solved and made on the
    GPU.  The energies of every row and the bytes of every dump must be those of ngpusph=0, on one rank and on
    several.  ngpusphcheck=1 makes the CPUs do the work too and stops the run at the first difference, so a run
    with it must also finish.  The CPU build has no GPU code: it must say so in log0.sph and give the same bits.
    """
    out = []
    for tag, exe in (('gpu', ctx.exe), ('cpu', ctx.cpu_exe)):
        if not exe:
            continue
        for np_ in (1, ctx.np):
            res = {}
            for mode, extra in (('0', dict(ngpusph=0)), ('1', dict(ngpusph=1)), ('1c', dict(ngpusph=1, ngpusphcheck=1))):
                d = ctx.dir('gs_%s_%d_%s' % (tag, np_, mode))
                p, e = run(exe, d, np=np_, **dict(POLY, tf=1, dtout=1, **extra))
                names = sorted(f for f in os.listdir(d) if f.startswith('out') and f.endswith('.sph'))
                if not names:
                    raise Fail('%s build, np=%d, ngpusph mode %s: no dump was written\n%s' % (tag, np_, mode, (p.stdout + p.stderr)[-800:]))
                log = open(os.path.join(d, 'log0.sph')).read() if os.path.exists(os.path.join(d, 'log0.sph')) else ''
                if mode != '0':
                    want = 'ngpusph: h, rho and neighbour lists on the CPUs: this build has no GPU code' if tag == 'cpu' else 'ngpusph: h, rho and neighbour lists on the GPU'
                    if want not in log:
                        raise Fail('%s build, np=%d, ngpusph mode %s: log0.sph does not say "%s"' % (tag, np_, mode, want))
                res[mode] = (energies(e), [(f, open(os.path.join(d, f), 'rb').read()) for f in names])
            for mode in ('1', '1c'):
                if res['0'][0] != res[mode][0]:
                    raise Fail('%s build, np=%d: ngpusph mode %s changes the energies' % (tag, np_, mode))
                if res['0'][1] != res[mode][1]:
                    raise Fail('%s build, np=%d: ngpusph mode %s changes a dump' % (tag, np_, mode))
            out.append('%s np=%d: %d rows, %d dumps' % (tag, np_, len(res['0'][0]), len(res['0'][1])))
    return 'ngpusph=0, 1 and 1 with the check give the same bits: ' + '; '.join(out)

def t_gpu_forces(ctx):
    """ngpusph=3 (the pair sums of the forces on the GPU too) agrees with ngpusph=1, to the rounding of the sums.

    Bit 2 of ngpusph puts div v, curl v, the accelerations and du/dt on the GPU, with the terms and the order of
    the sums of the CPU; they can differ from the CPU's in the last bit of a term (sqrt against pow in the
    viscosity).  An unrelaxed polytrope, whose particles approach each other and several time-step bins, is run
    with shared and with block steps: the energies of every row must agree with ngpusph=1 to one unit of the last
    digit that energy0.sph prints (7 digits: a rounding of a sum can change that digit, and 1e-9 asked for the
    same digits), and a run
    with ngpusphcheck=2 (the CPUs make the sums too and the differences go to log0.sph) must report
    differences of the accelerations below 1e-12 of the particle's and of du/dt below 1e-10 of the largest.
    The CPU build says that it has no GPU code and gives the same bits.
    """
    out = []
    for tag, exe in (('gpu', ctx.exe), ('cpu', ctx.cpu_exe)):
        if not exe:
            continue
        for np_ in (1, ctx.np):
            for blk in (0, 1):
                res, txt = {}, {}
                for mode, extra in (('1', dict(ngpusph=1)), ('3', dict(ngpusph=3)), ('3c', dict(ngpusph=3, ngpusphcheck=2))):
                    d = ctx.dir('gf_%s_%d_%d_%s' % (tag, np_, blk, mode))
                    p, e = run(exe, d, np=np_, **dict(POLY, nrelax=0, nblock=blk, tf=1, dtout=0.5, **extra))
                    log = open(os.path.join(d, 'log0.sph')).read() if os.path.exists(os.path.join(d, 'log0.sph')) else ''
                    if mode != '1':
                        want = ('ngpusph: the pair sums (forces) on the CPUs: this build has no GPU code' if tag == 'cpu'
                                else 'ngpusph: the pair sums (forces) on the GPU')
                        if want not in log:
                            raise Fail('%s build, np=%d, nblock=%d, ngpusph mode %s: log0.sph does not say "%s"\n%s'
                                       % (tag, np_, blk, mode, want, (p.stdout + p.stderr)[-600:]))
                    res[mode] = energies(e)
                    txt[mode] = energy_words(e)
                    if mode == '3c' and tag == 'gpu':
                        lines = [l for l in log.splitlines() if l.startswith(' ngpusphcheck forces:') or l.startswith('ngpusphcheck forces:')]
                        if not lines:
                            raise Fail('np=%d, nblock=%d: ngpusphcheck=2 wrote no line of differences' % (np_, blk))
                        worst = 0.0
                        worstdu = 0.0
                        for l in lines:
                            f = l.split()
                            a = float(f[f.index('a_rel_max') + 1])
                            du = float(f[f.index('du_max') + 1]) / max(float(f[f.index('du_ref') + 1]), 1e-300)
                            worst, worstdu = max(worst, a), max(worstdu, du)
                        if worst > 1e-12 or worstdu > 1e-10:
                            raise Fail('np=%d, nblock=%d: the forces on the GPU differ from the CPU by %.2e of a (limit 1e-12) '
                                       'and %.2e of du/dt (limit 1e-10)' % (np_, blk, worst, worstdu))
                for mode in ('3', '3c'):
                    if len(res['1']) != len(res[mode]):
                        raise Fail('%s build, np=%d, nblock=%d: ngpusph mode %s gives another number of rows' % (tag, np_, blk, mode))
                    for r1, r3, w1, w3 in zip(res['1'], res[mode], txt['1'], txt[mode]):
                        for k, name in ((1, 'W'), (3, 'U'), (4, 'Etot')):
                            if tag == 'cpu':
                                if r1[k] != r3[k]:
                                    raise Fail('%s build, np=%d, nblock=%d: ngpusph mode %s changes %s' % (tag, np_, blk, mode, name))
                            else:
                                close_digits(w3[k], w1[k], 1, '%s, np=%d, nblock=%d, ngpusph=%s against 1' % (name, np_, blk, mode))
                out.append('%s np=%d nblock=%d: %d rows' % (tag, np_, blk, len(res['1'])))
    return 'ngpusph=3 agrees with 1 (and the check of the forces is clean): ' + '; '.join(out)

def t_gpu_ranks(ctx):
    """ngpusphcheck=2 and 3 with more ranks than GPUs: a rank without a GPU of its own takes part in the comparison (it left the run hanging).

    With ngpuperrank=1 rank r uses GPU r, and a rank beyond the GPUs of the node does the SPH on the CPUs.  The comparison of the forces
    (ngpusphcheck=2) adds the sums of all the ranks with collective calls, which such a rank used to skip.  CUDA_VISIBLE_DEVICES=0 shows one GPU only, so
    with 2 or 3 ranks one or two of them have none, on any node.  An unrelaxed polytrope is run with ngpusph=3 and ngpusphcheck=2 and 3 (and with
    ngpusph=1): it must finish within 5 minutes, write the lines of the comparison and give the energies of ngpusph=1 to 1e-9 (the check goes on with the
    forces of the CPUs).  The CPU build has no GPUs and is not run.
    """
    if not ctx.exe:
        return 'no GPU executable'
    out = []
    old = os.environ.get('CUDA_VISIBLE_DEVICES')
    os.environ['CUDA_VISIBLE_DEVICES'] = '0'
    try:
        for np_ in (2, 3):
            res = {}
            for mode, extra in (('1', dict(ngpusph=1)), ('3c', dict(ngpusph=3, ngpusphcheck=2)), ('3cc', dict(ngpusph=3, ngpusphcheck=3))):
                d = ctx.dir('gk_%d_%s' % (np_, mode))
                try:
                    p, e = run(ctx.exe, d, np=np_, timeout=300, **dict(POLY, nrelax=0, nblock=0, tf=1, dtout=0.5, **extra))
                except subprocess.TimeoutExpired:
                    raise Fail('np=%d, ngpusph mode %s: no end after 5 minutes (a rank waits for the others)' % (np_, mode))
                res[mode] = energies(e)
                if mode != '1':
                    log = open(os.path.join(d, 'log0.sph')).read()
                    lines = [l for l in log.splitlines() if 'ngpusphcheck forces:' in l]
                    if not lines:
                        raise Fail('np=%d, ngpusph mode %s: no line of differences in log0.sph\n%s' % (np_, mode, (p.stdout + p.stderr)[-600:]))
                    for l in lines:
                        f = l.split()
                        if float(f[f.index('a_rel_max') + 1]) > 1e-12:
                            raise Fail('np=%d, ngpusph mode %s: the forces on the GPU differ from the CPU: %s' % (np_, mode, l.strip()))
            for mode in ('3c', '3cc'):
                if len(res['1']) != len(res[mode]):
                    raise Fail('np=%d, ngpusph mode %s gives another number of rows' % (np_, mode))
                for r1, r3 in zip(res['1'], res[mode]):
                    for k, name in ((1, 'W'), (3, 'U'), (4, 'Etot')):
                        close(r3[k], r1[k], 1e-9, '%s, np=%d, ngpusph mode %s against 1' % (name, np_, mode))
            out.append('np=%d: %d rows' % (np_, len(res['1'])))
    finally:
        if old is None:
            os.environ.pop('CUDA_VISIBLE_DEVICES', None)
        else:
            os.environ['CUDA_VISIBLE_DEVICES'] = old
    return 'ngpusphcheck with ranks that have no GPU of their own finishes and agrees with ngpusph=1: ' + '; '.join(out)

def t_gpu_resident(ctx):
    """ngpusph=7 (the candidates of the smoothing-length search kept on the GPU between steps: Verlet lists) gives exactly the bits of ngpusph=3.

    An unrelaxed polytrope, whose particles move, is run with the margin of the lists made small (ngpusphskin=0.05) so that they are built again several times in the run.
    The energies of every row and the bytes of every dump must be those of ngpusph=3, on one rank and on several; and with ngpusphcheck=3 (the GPU against the CPU at every
    step: the smoothing lengths, the lists and the forces, going on with the CPU forces) those of ngpusph=3 with ngpusphcheck=3.  The CPU build says that it has no GPU code and gives the same bits.
    """
    out = []
    for tag, exe in (('gpu', ctx.exe), ('cpu', ctx.cpu_exe)):
        if not exe:
            continue
        for np_ in (1, ctx.np):
            res = {}
            for mode, extra in (('3', dict(ngpusph=3)), ('7', dict(ngpusph=7, ngpusphskin='0.05d0')), ('3c', dict(ngpusph=3, ngpusphcheck=3)),
                                ('7c', dict(ngpusph=7, ngpusphskin='0.05d0', ngpusphcheck=3))):
                d = ctx.dir('gr_%s_%d_%s' % (tag, np_, mode))
                p, e = run(exe, d, np=np_, **dict(POLY, nrelax=0, tf=1, dtout=1, **extra))
                names = sorted(f for f in os.listdir(d) if f.startswith('out') and f.endswith('.sph'))
                if not names:
                    raise Fail('%s build, np=%d, ngpusph mode %s: no dump was written\n%s' % (tag, np_, mode, (p.stdout + p.stderr)[-800:]))
                log = open(os.path.join(d, 'log0.sph')).read() if os.path.exists(os.path.join(d, 'log0.sph')) else ''
                want = 'ngpusph: h, rho and neighbour lists on the CPUs: this build has no GPU code' if tag == 'cpu' else 'ngpusph: h, rho and neighbour lists on the GPU'
                if want not in log:
                    raise Fail('%s build, np=%d, ngpusph mode %s: log0.sph does not say "%s"' % (tag, np_, mode, want))
                res[mode] = (energies(e), [(f, open(os.path.join(d, f), 'rb').read()) for f in names])
            # (with ngpusphcheck the run goes on with the forces of the CPU, so the runs with the check are compared with each other)
            for mode, ref in (('7', '3'), ('7c', '3c')):
                if res[ref][0] != res[mode][0]:
                    raise Fail('%s build, np=%d: ngpusph mode %s changes the energies of mode %s' % (tag, np_, mode, ref))
                if res[ref][1] != res[mode][1]:
                    raise Fail('%s build, np=%d: ngpusph mode %s changes a dump of mode %s' % (tag, np_, mode, ref))
            out.append('%s np=%d: %d rows, %d dumps' % (tag, np_, len(res['3'][0]), len(res['3'][1])))
    return 'ngpusph=7 gives the bits of 3, also with the check against the CPU: ' + '; '.join(out)

def t_gpu_mixed(ctx):
    """ngpusphprec=1 (the terms of the pairs of the forces on the GPU in single precision, sums in double) stays within the tolerances given below.

    An unrelaxed polytrope is run with ngpusph=3, ngpusphprec=1 and ngpusphcheck=2 (the CPUs make the forces too and the differences go to log0.sph; the run goes on with the CPU
    forces): the accelerations must differ by less than 1e-4 of the particle's and du/dt by less than 1e-3 of the largest, at every step.  The same run without the check must end with
    energies within 1e-4 of those of ngpusph=3 in double precision.  The CPU build says it has no GPU code and gives the same bits.
    """
    out = []
    for tag, exe in (('gpu', ctx.exe), ('cpu', ctx.cpu_exe)):
        if not exe:
            continue
        for np_ in (1, ctx.np):
            res = {}
            for mode, extra in (('d', dict(ngpusph=3)), ('m', dict(ngpusph=3, ngpusphprec=1)), ('mc', dict(ngpusph=3, ngpusphprec=1, ngpusphcheck=2))):
                d = ctx.dir('gm_%s_%d_%s' % (tag, np_, mode))
                p, e = run(exe, d, np=np_, **dict(POLY, nrelax=0, tf=1, dtout=1, **extra))
                log = open(os.path.join(d, 'log0.sph')).read() if os.path.exists(os.path.join(d, 'log0.sph')) else ''
                res[mode] = energies(e)
                if mode == 'mc' and tag == 'gpu':
                    lines = [l for l in log.splitlines() if 'ngpusphcheck forces:' in l]
                    if not lines:
                        raise Fail('np=%d: ngpusphcheck=2 wrote no line of differences' % np_)
                    worst = worstdu = 0.0
                    for l in lines:
                        f = l.split('forces:')[1].split()
                        worst = max(worst, float(f[f.index('a_rel_max') + 1]))
                        worstdu = max(worstdu, float(f[f.index('du_max') + 1]) / max(float(f[f.index('du_ref') + 1]), 1e-300))
                    if worst > 1e-4 or worstdu > 1e-3:
                        raise Fail('np=%d: single precision differs from the CPU by %.2e of a (limit 1e-4) and %.2e of du/dt (limit 1e-3)' % (np_, worst, worstdu))
                    out.append('np=%d: a %.1e, du/dt %.1e' % (np_, worst, worstdu))
            tol = 0.0 if tag == 'cpu' else 1e-4
            if len(res['d']) != len(res['m']):
                raise Fail('%s build, np=%d: another number of rows' % (tag, np_))
            for r1, r2 in zip(res['d'], res['m']):
                for k, name in ((1, 'W'), (3, 'U'), (4, 'Etot')):
                    if tol == 0.0:
                        if r1[k] != r2[k]:
                            raise Fail('%s build, np=%d: ngpusphprec=1 changes %s' % (tag, np_, name))
                    else:
                        close(r2[k], r1[k], tol, '%s, np=%d, ngpusphprec=1 against 0' % (name, np_))
    return 'ngpusphprec=1 within the tolerances' + (': ' + '; '.join(out) if out else ' (CPU build: the same bits)')

def t_sph_auto(ctx):
    """ngpusph=-1 chooses from the configuration, says so in log0.sph, and gives exactly the bits of the value it chose, set by hand.

    The line "SPH on the GPU: ngpusph=-1 (auto): <why> -> ngpusph=N" must be in log0.sph; the energies of every row and the bytes of every dump must be those of a run with ngpusph=N.  The CPU build must
    choose 0.  With a GPU build and a tree on the GPU or on the CPUs the choice is 7; with the direct sum it depends on the threads per rank, the particles (a star below and one above each limit) and the ranks
    (several ranks: 1).
    """
    out = []
    for tag, exe in (('gpu', ctx.exe), ('cpu', ctx.cpu_exe)):
        if not exe:
            continue
        cpu = exe.endswith('_cpu_sph')
        if cpu:
            for np_ in (1, ctx.np):
                for gname, gx in (('direct', dict(ngravtree=0)), ('tree', dict(ngravtree=1, treeacc='0.5d0'))):
                    chosen = auto_case(ctx, exe, 'sa_cpu_%d_%s' % (np_, gname), np_, 1, 5000, gx, (None, 0))
                    out.append('cpu np=%d %s -> %d' % (np_, gname, chosen[1]))
            continue
        for label, np_, threads, nn, written, want in AUTO_CASES_SPH:
            chosen = auto_case(ctx, exe, 'sa_' + label, np_, threads, nn, written, (None, want))
            if chosen is None:
                out.append('%s skipped (no OpenMP)' % label)
                continue
            out.append('%s -> %d' % (label, chosen[1]))
    return 'ngpusph=-1 chooses as the rule says, says so, and gives the bits of its choice: ' + '; '.join(out)

def t_vdotsm(ctx):
    """With cn2 and cn4 finite (tstep uses the smoothed accelerations of vdotsm) ngpusph=1 gives the bits of 0 and ngpusph=7 those of 3.

    vdotsm sums the accelerations over the neighbour lists in list(), so whoever keeps the lists on the GPU (ngpusph=7) has to copy them when the courant numbers cn2 or cn4
    make tstep read the result; before, vdotsm read lists of an older step and the time step, and the run, came out wrong.  An unrelaxed polytrope (nrelax=0) is run with
    cn2=0.06 and cn4=0.1 and ngpusph 0, 1, 3 and 7: the energies of every row, the bytes of every dump and the number of steps of ngpusph=1 must be those of 0, and of 7 those of 3
    (the CPU build: all four the same).  The same run with cn2 and cn4 at their default (1d30) must come out different, or the test would prove nothing.
    """
    out = []
    for tag, exe in (('gpu', ctx.exe), ('cpu', ctx.cpu_exe)):
        if not exe:
            continue
        for np_ in (1, ctx.np):
            res = {}
            for mode, extra in (('0', dict(ngpusph=0)), ('1', dict(ngpusph=1)), ('3', dict(ngpusph=3)), ('7', dict(ngpusph=7)),
                                ('d', dict(ngpusph=7, cn2='1.d30', cn4='1.d30'))):
                d = ctx.dir('vd_%s_%d_%s' % (tag, np_, mode))
                p, e = run(exe, d, np=np_, **dict(POLY, nrelax=0, tf=1, dtout=1, **dict(dict(cn2='0.06d0', cn4='0.1d0'), **extra)))
                names = sorted(f for f in os.listdir(d) if f.startswith('out') and f.endswith('.sph'))
                if not names:
                    raise Fail('%s build, np=%d, mode %s: no dump was written\n%s' % (tag, np_, mode, (p.stdout + p.stderr)[-800:]))
                steps = len(re.findall(r'^mainit', p.stdout, re.M))
                res[mode] = (energies(e), [(f, open(os.path.join(d, f), 'rb').read()) for f in names], steps)
            for mode, ref in (('1', '0'), ('7', '3')) + ((('3', '0'), ('7', '0')) if tag == 'cpu' else ()):
                for k, what in ((0, 'energies'), (1, 'a dump'), (2, 'the number of steps')):
                    if res[mode][k] != res[ref][k]:
                        raise Fail('%s build, np=%d: ngpusph=%s changes %s of ngpusph=%s (cn2=0.06, cn4=0.1, nrelax=0)' % (tag, np_, mode, what, ref))
            if res['7'][2] == res['d'][2] and res['7'][0] == res['d'][0]:
                raise Fail('%s build, np=%d: cn2=0.06 and cn4=0.1 do not change the run (%d steps): the test proves nothing' % (tag, np_, res['7'][2]))
            out.append('%s np=%d: %d steps (%d with the default cn2, cn4)' % (tag, np_, res['7'][2], res['d'][2]))
    return 'with cn2 and cn4 finite ngpusph=1 gives the bits of 0 and 7 those of 3: ' + '; '.join(out)

def t_gpu_multidev(ctx):
    """The SPH and the forces with several devices for the rank (SPHGPU_NDEV: logical devices dealt out round the GPUs) run to the end and agree with one device.

    SPHGPU_NDEV=K (a variable of the tests, see docs/source/using/gpu_sph.rst) makes the SPH use K logical devices on the GPU or GPUs of the rank, each with its own stream and buffers: one
    GPU then gives the code paths and the sums of K devices (the times are not those of K GPUs).  The unrelaxed polytrope is run with ngpusph=3 and ngpusphcheck=2 with 1, 2 and 4 logical
    devices: every run must finish, log0.sph must say that the forces are on the GPU, SPHGPU_PROF=1 must say that the devices were 1, 2 and 4, the differences that ngpusphcheck=2 writes
    must be below 1e-12 of the acceleration of a particle and 1e-10 of du/dt, and the energies of every row must be those of one device to 1e-9.  With 3 logical devices and ngpusph=7 the
    dump and the energies must not depend on SPHGPU_FLAGMOD (7, 2 or none): the particles that the CPU solves in place of the GPU do not change how the forces are summed.
    """
    out = []
    if not ctx.exe:
        raise Fail('no GPU executable')
    res = {}
    saved = {k: os.environ.get(k) for k in ('SPHGPU_NDEV', 'SPHGPU_PROF')}
    try:
        for nd in (1, 2, 4):
            os.environ['SPHGPU_NDEV'] = str(nd)
            os.environ['SPHGPU_PROF'] = '1'
            d = ctx.dir('md_%d' % nd)
            p, e = run(ctx.exe, d, np=1, **dict(POLY, nrelax=0, tf=1, dtout=0.5, ngpusph=3, ngpusphcheck=2))
            log = open(os.path.join(d, 'log0.sph')).read() if os.path.exists(os.path.join(d, 'log0.sph')) else ''
            if 'ngpusph: the pair sums (forces) on the GPU' not in log:
                raise Fail('%d logical devices: log0.sph does not say that the forces are on the GPU\n%s' % (nd, (p.stdout + p.stderr)[-600:]))
            if not re.search(r'%d devices' % nd, p.stderr):
                raise Fail('%d logical devices: SPHGPU_PROF=1 does not report %d devices' % (nd, nd))
            lines = [l for l in log.splitlines() if 'ngpusphcheck forces:' in l]
            if not lines:
                raise Fail('%d logical devices: ngpusphcheck=2 wrote no line of differences' % nd)
            worst = worstdu = 0.0
            for l in lines:
                f = l.split()
                a = float(f[f.index('a_rel_max') + 1])
                du = float(f[f.index('du_max') + 1]) / max(float(f[f.index('du_ref') + 1]), 1e-300)
                worst, worstdu = max(worst, a), max(worstdu, du)
            if worst > 1e-12 or worstdu > 1e-10:
                raise Fail('%d logical devices: the forces on the GPU differ from the CPU by %.2e of a (limit 1e-12) and %.2e of du/dt (limit 1e-10)' % (nd, worst, worstdu))
            res[nd] = energies(e)
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    for nd in (2, 4):
        if len(res[nd]) != len(res[1]):
            raise Fail('%d logical devices give another number of rows than one' % nd)
        for r1, rk in zip(res[1], res[nd]):
            for k, name in ((1, 'W'), (3, 'U'), (4, 'Etot')):
                close(rk[k], r1[k], 1e-9, '%s with %d logical devices against 1' % (name, nd))
        out.append('%d devices: %d rows' % (nd, len(res[nd])))
    # With several devices a step with particles the CPU solves (more than 128 candidates; SPHGPU_FLAGMOD=K leaves every Kth particle to the CPU as if they had) is summed as the
    # other steps are: the results do not depend on which steps have such particles.  ngpusph=7, no check against the CPU: the bytes of the dump and the energies of every row.
    ref = None
    saved = {k: os.environ.get(k) for k in ('SPHGPU_NDEV', 'SPHGPU_FLAGMOD')}
    try:
        os.environ['SPHGPU_NDEV'] = '3'
        for fm in (0, 7, 2):
            os.environ['SPHGPU_FLAGMOD'] = str(fm)
            d = ctx.dir('md_fm%d' % fm)
            p, e = run(ctx.exe, d, np=1, **dict(POLY, nrelax=0, tf=1, dtout=1, ngpusph=7))
            names = sorted(f for f in os.listdir(d) if f.startswith('out') and f.endswith('.sph'))
            if not names:
                raise Fail('3 logical devices, SPHGPU_FLAGMOD=%d: no dump was written\n%s' % (fm, (p.stdout + p.stderr)[-800:]))
            got = (energies(e), [(f, open(os.path.join(d, f), 'rb').read()) for f in names])
            if ref is None:
                ref = got
            elif got != ref:
                raise Fail('3 logical devices: SPHGPU_FLAGMOD=%d changes the energies or a dump: the sums depend on which particles the CPU solved' % fm)
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    out.append('SPHGPU_FLAGMOD 7 and 2 give the bits of none')
    return 'logical devices run to the end, the check of the forces is clean, the energies are those of one device, and the particles solved by the CPU do not change the sums: ' + '; '.join(out)

# --------------------------------------------------------------- HDF5 -----

def t_hdf5_dump(ctx):
    """nhdf5=1 writes an out*.h5 that holds exactly what out*.sph holds.

    The two files are compared bit for bit, column by column (the floats
    through their bit patterns) and entry by entry of the header.  nhdf5=2 must
    write the same out*.h5 and no out*.sph, and restartrad.sph stays native.
    A build without the HDF5 writer stops a run that asks for nhdf5>0 (in init,
    with a message), and the test is skipped then, as it is when h5py is
    missing.  Build with make HDF5=1, or run the suite with --hdf5, to run it.
    """
    import glob, gzip
    try:
        import numpy, h5py                                   # noqa: F401
    except ImportError as e:
        raise Skip('needs numpy and h5py (%s)' % e)
    sys.path.insert(0, os.path.join(ROOT, 'parallel_bleeding_edge', 'bin'))
    import StSm_snapshot as ss

    s = dict(POLY, tf=2, dtout=1)
    d = ctx.dir('hdf5_dump')
    p, _ = run(ctx.exe, d, np=ctx.np, nhdf5=1, **s)
    text = p.stdout + p.stderr
    if 'needs a build with make HDF5=1' in text:
        raise Skip('this build has no HDF5 writer; build with make HDF5=1 '
                   '(tests/run_tests.py --hdf5)')
    sph = sorted(glob.glob(os.path.join(d, 'out*.sph')))
    h5 = sorted(glob.glob(os.path.join(d, 'out*.h5')))
    if len(sph) < 2:
        raise Fail('the run wrote %d snapshots, expected at least 2\n' % len(sph)
                   + text[-800:])
    if [f[:-4] for f in sph] != [f[:-3] for f in h5]:
        raise Fail('nhdf5=1 wrote %s but %s' % ([os.path.basename(f) for f in h5],
                                                 [os.path.basename(f) for f in sph]))
    for a, b in zip(sph, h5):
        nd, bad, badh = ss.compare(a, b)
        if nd != len(ss.COLUMNS) or bad or badh:
            raise Fail('%s and %s differ: %d columns compared, columns %s, '
                       'header entries %s' % (os.path.basename(a),
                                              os.path.basename(b), nd, bad, badh))
    if os.path.exists(os.path.join(d, 'restartrad.sph')):
        if ss.is_hdf5(os.path.join(d, 'restartrad.sph')):
            raise Fail('restartrad.sph is not in the native format')
    log = open(os.path.join(d, 'log0.sph')).read()
    if len(re.findall(r'wall time writing out\d+\.h5 \(hdf5\)', log)) != len(h5):
        raise Fail('log0.sph does not give the time of each HDF5 write')

    # nhdf5=2: the same .h5 files, and no .sph
    d2 = ctx.dir('hdf5_dump_only')
    run(ctx.exe, d2, np=ctx.np, nhdf5=2, **s)
    if glob.glob(os.path.join(d2, 'out*.sph')):
        raise Fail('nhdf5=2 still wrote out*.sph')
    h5b = sorted(glob.glob(os.path.join(d2, 'out*.h5')))
    if [os.path.basename(f) for f in h5b] != [os.path.basename(f) for f in h5]:
        raise Fail('nhdf5=2 wrote %s, nhdf5=1 wrote %s'
                   % ([os.path.basename(f) for f in h5b],
                      [os.path.basename(f) for f in h5]))
    for a, b in zip(sph, h5b):
        nd, bad, badh = ss.compare(a, b)
        if bad or badh:
            raise Fail('nhdf5=2 gave %s, which differs from the nhdf5=1 run\'s '
                       '%s in columns %s, header entries %s'
                       % (os.path.basename(b), os.path.basename(a), bad, badh))

    # the other layout, with ueq and tthermal, needs the cooling tables
    gz = os.path.join(ROOT, 'example_input', 'sph.opacity.gz')
    extra = ''
    if os.path.exists(gz):
        d3 = ctx.dir('hdf5_dump_cooling')
        with gzip.open(gz, 'rb') as fi, open(os.path.join(d3, 'sph.opacity'), 'wb') as fo:
            shutil.copyfileobj(fi, fo)
        run(ctx.exe, d3, np=ctx.np, nhdf5=1, ncooling=1, **s)
        sph3 = sorted(glob.glob(os.path.join(d3, 'out*.sph')))
        h53 = sorted(glob.glob(os.path.join(d3, 'out*.h5')))
        if not sph3 or len(sph3) != len(h53):
            raise Fail('the ncooling=1 run wrote %d .sph and %d .h5'
                       % (len(sph3), len(h53)))
        for a, b in zip(sph3, h53):
            nd, bad, badh = ss.compare(a, b)
            if nd != len(ss.COLUMNS) + 2 or bad or badh:
                raise Fail('ncooling=1: %s differs from %s (%d columns, %s, %s)'
                           % (os.path.basename(b), os.path.basename(a), nd, bad, badh))
        extra = ', and with ncooling=1 (20 columns)'
    return ('%d snapshots: out*.h5 equals out*.sph bit for bit (18 columns, 25 '
            'header entries); nhdf5=2 gives the same files without out*.sph%s'
            % (len(sph), extra))

TESTS = {
    'virial':           t_virial,
    'energy':           t_energy,
    'quiet':            t_quiet,
    'ranks':            t_ranks,
    'cpu_gpu':          t_cpu_gpu,
    'block_cpu_gpu':    t_block_cpu_gpu,
    'gpu_hilo':         t_gpu_hilo,
    'gpu_tree':         t_gpu_tree,
    'gravity_auto':     t_gravity_auto,
    'gpu_sph':          t_gpu_sph,
    'gpu_forces':       t_gpu_forces,
    'gpu_ranks':        t_gpu_ranks,
    'gpu_resident':     t_gpu_resident,
    'gpu_mixed':        t_gpu_mixed,
    'sph_auto':         t_sph_auto,
    'vdotsm':           t_vdotsm,
    'gpu_multidev':     t_gpu_multidev,
    'guard_ngravprocs': t_guard_ngravprocs,
    'startfiles':       t_startfiles,
    'periodic':         t_periodic,
    'periodic_block':   t_periodic_block,
    'guard_periodic':   t_guard_periodic,
    'hdf5_dump':        t_hdf5_dump,
}

# ------------------------------------------------------------------- main ---

class Ctx:
    def __init__(self, base, exe, cpu_exe, np):
        self.base, self.exe, self.cpu_exe, self.np = base, exe, cpu_exe, np
    def dir(self, name):
        d = os.path.join(self.base, name)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        return d

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('names', nargs='*', help='tests to run (default: all)')
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--np', type=int, default=4, help='MPI ranks (default 4)')
    ap.add_argument('--exe', help='use this executable instead of building')
    ap.add_argument('--cpu-exe', help='CPU executable for the cpu_gpu test')
    ap.add_argument('--with-cpu', action='store_true', help='also build the CPU version')
    ap.add_argument('--hdf5', action='store_true',
                    help='build with make HDF5=1, which writes out*.h5 (the hdf5_dump test)')
    ap.add_argument('--keep', action='store_true', help='keep the scratch directory')
    ap.add_argument('--quick', action='store_true',
                    help='use a smaller star, for continuous integration')
    a = ap.parse_args()

    # The assertions are unchanged by --quick; only the particle count is.
    # A 2000-particle polytrope reproduces the virial theorem to the same
    # four figures as a 10000-particle one and takes about a fifteenth of
    # the time, which is the difference between a usable CI job and one
    # nobody waits for.
    if a.quick:
        POLY['n'] = 2000

    if a.list:
        for n, f in TESTS.items():
            print('%-18s %s' % (n, (f.__doc__ or '').strip().split('\n')[0]))
        return 0

    names = a.names or list(TESTS)
    unknown = [n for n in names if n not in TESTS]
    if unknown:
        print('unknown test(s): %s' % ', '.join(unknown)); return 2

    base = tempfile.mkdtemp(prefix='starsmasher-tests-')
    print('scratch: %s' % base)
    try:
        exe = a.exe or build('gpu', base, a.hdf5)
        cpu = a.cpu_exe or (build('cpu', base, a.hdf5) if a.with_cpu else None)
        ctx = Ctx(base, exe, cpu, a.np)
        print('executable: %s' % exe)
        print('ranks: %d\n' % a.np)

        width = max(len(n) for n in names)
        failures = skipped = 0
        for n in names:
            sys.stdout.write('%-*s ... ' % (width, n)); sys.stdout.flush()
            t0 = time.time()
            try:
                detail = TESTS[n](ctx)
                print('ok   (%4.1fs)  %s' % (time.time() - t0, detail))
            except Skip as exc:
                skipped += 1
                print('skip (%4.1fs)  %s' % (time.time() - t0, exc))
            except Fail as exc:
                failures += 1
                print('FAIL (%4.1fs)' % (time.time() - t0))
                for line in str(exc).split('\n'):
                    print('    %s' % line)
            except subprocess.TimeoutExpired:
                failures += 1
                print('TIMEOUT')
        print('\n%d of %d passed%s' % (len(names) - failures - skipped, len(names),
              ', %d skipped' % skipped if skipped else ''))
        return 1 if failures else 0
    finally:
        if not a.keep:
            shutil.rmtree(base, ignore_errors=True)
        else:
            print('kept: %s' % base)

if __name__ == '__main__':
    sys.exit(main())
