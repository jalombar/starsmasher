#!/usr/bin/env python3
"""StarSmasher test suite.

Runs short calculations and checks them against expectations that must hold,
rather than against stored output.  See tests/README.md.
"""
import argparse, math, os, re, shutil, subprocess, sys, tempfile, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(ROOT, 'parallel_bleeding_edge', 'src')

class Fail(Exception):
    pass

# ---------------------------------------------------------------- helpers ---

def sh(cmd, cwd=None, timeout=1800, check=True):
    p = subprocess.run(cmd, cwd=cwd, shell=isinstance(cmd, str),
                       capture_output=True, text=True, timeout=timeout)
    if check and p.returncode != 0:
        raise Fail('command failed: %s\n%s' % (cmd, (p.stderr or p.stdout)[-1500:]))
    return p

def build(target='gpu', base=None):
    """Build in a copy of the source tree, never in the repository itself.

    Building in place would leave object files, .mod files and executables
    scattered through parallel_bleeding_edge/, which is both untidy and a
    trap for anyone who then runs `git add -A`.
    """
    work = os.path.join(base, 'build-%s' % target, 'src')
    shutil.copytree(SRC, work)
    p = sh('make %s' % ('' if target == 'gpu' else target), cwd=work, check=False)
    parent = os.path.dirname(work)
    exe = [os.path.join(parent, f) for f in os.listdir(parent)
           if f.endswith('_%s_sph' % target)]
    if not exe:
        raise Fail('build produced no _%s_sph executable\n%s'
                   % (target, (p.stdout + p.stderr)[-1500:]))
    return exe[0]

def sph_input(**kw):
    body = '\n'.join(' %s=%s,' % (k, v) for k, v in kw.items())
    return ' &input\n%s\n &end\n' % body

def run(exe, workdir, np=4, init="'1es'", timeout=1800, **settings):
    os.makedirs(workdir, exist_ok=True)
    shutil.copy(exe, os.path.join(workdir, 'run_sph'))
    open(os.path.join(workdir, 'sph.init'), 'w').write(
        ' &INITT\n INAME=%s\n &END\n' % init)
    open(os.path.join(workdir, 'sph.input'), 'w').write(sph_input(**settings))
    p = sh(['mpirun', '--oversubscribe', '-np', str(np), './run_sph'],
           cwd=workdir, timeout=timeout, check=False)
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

def t_cpu_gpu(ctx):
    """The CPU and GPU builds must agree."""
    if not ctx.cpu_exe:
        raise Fail('no CPU build available (pass --with-cpu to build one)')
    res = {}
    for tag, exe in (('gpu', ctx.exe), ('cpu', ctx.cpu_exe)):
        _, e = run(exe, ctx.dir('cpugpu_' + tag), np=ctx.np, **POLY)
        res[tag] = energies(e)[0]
    for i, name in ((1, 'W'), (3, 'U'), (4, 'Etot')):
        close(res['cpu'][i], res['gpu'][i], 1e-6, '%s, cpu vs gpu' % name)
    return 'cpu and gpu agree on W, U and Etot to 1e-6'

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

TESTS = {
    'virial':           t_virial,
    'energy':           t_energy,
    'quiet':            t_quiet,
    'ranks':            t_ranks,
    'cpu_gpu':          t_cpu_gpu,
    'guard_ngravprocs': t_guard_ngravprocs,
    'startfiles':       t_startfiles,
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
        exe = a.exe or build('gpu', base)
        cpu = a.cpu_exe or (build('cpu', base) if a.with_cpu else None)
        ctx = Ctx(base, exe, cpu, a.np)
        print('executable: %s' % exe)
        print('ranks: %d\n' % a.np)

        width = max(len(n) for n in names)
        failures = 0
        for n in names:
            sys.stdout.write('%-*s ... ' % (width, n)); sys.stdout.flush()
            t0 = time.time()
            try:
                detail = TESTS[n](ctx)
                print('ok   (%4.1fs)  %s' % (time.time() - t0, detail))
            except Fail as exc:
                failures += 1
                print('FAIL (%4.1fs)' % (time.time() - t0))
                for line in str(exc).split('\n'):
                    print('    %s' % line)
            except subprocess.TimeoutExpired:
                failures += 1
                print('TIMEOUT')
        print('\n%d of %d passed' % (len(names) - failures, len(names)))
        return 1 if failures else 0
    finally:
        if not a.keep:
            shutil.rmtree(base, ignore_errors=True)
        else:
            print('kept: %s' % base)

if __name__ == '__main__':
    sys.exit(main())
