#!/usr/bin/env python3
"""Generate the sph.input reference page from parallel_bleeding_edge/src/init.f.

The namelist members, their default values and the explanatory comments all
come straight from the source, so this page cannot drift from the code.  Run
    python3 docs/generate_reference.py
from the top of the repository after changing init.f.
"""
import io, os, re, sys, textwrap

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INIT = os.path.join(ROOT, 'parallel_bleeding_edge', 'src', 'init.f')
OUT  = os.path.join(ROOT, 'docs', 'source', 'reference', 'sph_input.rst')

def rst_escape(text):
    """Comments in init.f contain characters RST treats as markup.

    ``|`` starts a substitution reference and ``*`` starts emphasis, so a
    comment such as ``dt2=cn2*(h/|a-a_smoothed|)^0.5`` would otherwise break
    the build.  Escape them rather than reformatting the comment, so the page
    shows exactly what the source says.
    """
    for ch in ('\\', '|', '*', '`', '_'):
        text = text.replace(ch, '\\' + ch)
    return text


def is_comment(line):
    """Fixed-form Fortran: a comment marker must sit in column 1."""
    return (not line.strip()) or line[0] in 'cC*!'


def applies_to(src_dir, members):
    """Work out, for each setting, which initialization scripts read it.

    A setting the integrator reads means the same thing for every INAME.  One
    that only a setup routine reads is specialised, and naming which routines
    saves a reader wondering why it appears to do nothing.

    The unit is the file, not the subroutine: initialize_parent.f holds helper
    routines besides `parent` itself, and a variable one of those reads is still
    part of building that star.  init.f is the exception, since it holds both
    the namelist plumbing and `readin`, so it is split by subroutine.
    """
    import glob

    def strip(t):
        return '\n'.join(l for l in t.split('\n') if l[:1] not in 'cC*!')

    init_path = os.path.join(src_dir, 'init.f')
    init_src = io.open(init_path, encoding='utf-8').read()
    disp = re.findall(r"iname\.eq\.'([a-z0-9]{3})'\s*\)\s*then\s*\n\s*call\s+(\w+)",
                      init_src)

    # which file defines each dispatched routine
    where, subtext = {}, {}
    for f in sorted(glob.glob(os.path.join(src_dir, '*.f')) +
                    glob.glob(os.path.join(src_dir, '*.f90'))):
        t = io.open(f, encoding='utf-8', errors='replace').read()
        parts = re.split(r'(?im)^[ \t]{0,8}(?:subroutine|program)\s+(\w+)', t)
        for i in range(1, len(parts), 2):
            n = parts[i].lower()
            where.setdefault(n, os.path.basename(f))
            subtext[n] = subtext.get(n, '') + strip(parts[i + 1])

    file_codes, sub_codes = {}, {}
    for code, sub in disp:
        sub = sub.lower()
        f = where.get(sub)
        if f == 'init.f':                 # rin lives beside the plumbing
            sub_codes.setdefault(sub, []).append(code)
        elif f:
            file_codes.setdefault(f, []).append(code)

    text = {}
    for f in sorted(glob.glob(os.path.join(src_dir, '*.f')) +
                    glob.glob(os.path.join(src_dir, '*.f90'))):
        text[os.path.basename(f)] = strip(
            io.open(f, encoding='utf-8', errors='replace').read())

    # Everything that is not a setup file, and not the namelist plumbing.
    evolution = '\n'.join(t for f, t in text.items()
                          if f not in file_codes and f != 'init.f')

    def seen(body, v):
        return re.search(r'(?<![A-Za-z0-9_])%s(?![A-Za-z0-9_])' % re.escape(v),
                         body, re.I) is not None

    out = {}
    for v in members:
        codes = []
        for f, cs in file_codes.items():
            if seen(text.get(f, ''), v):
                codes.extend(cs)
        for sub, cs in sub_codes.items():
            if seen(subtext.get(sub, ''), v):
                codes.extend(cs)
        out[v.lower()] = (sorted(set(codes)), seen(evolution, v))
    return out

SCRIPT_TITLES = {
    '1es': 'a single polytrope',
    '1mc': 'a single polytrope with a compact core',
    '2cr': 'a corotating binary',
    'bps': 'a binary plus a single star',
    'bph': 'a binary plus a black hole',
    'hyp': 'two bodies on a Keplerian orbit',
    'hbs': 'a binary encountering a single star',
    'erg': 'a star from a stellar-evolution profile',
    'meq': 'a star whose particle-mass scheme varies with position',
    'res': 'a rescaling of an existing model',
    'tri': 'a triple system',
    'bhe': 'an encounter with a supermassive black hole',
    'grs': 'a general-relativistic setup',
    'rin': 'a restart from an existing dump',
    'txt': 'particles laid out from an ASCII image',
}


def write_by_script(out_path, members, defaults, usage):
    """The same facts arranged the other way round: per script, what it reads.

    The sph.input page answers "where does this setting apply".  Someone who has
    already chosen a script has the opposite question, and answering it from a
    75-row table means reading every row.
    """
    o = []
    w = o.append
    w('.. This page is generated by docs/generate_reference.py.  Do not edit it.')
    w('')
    w('Settings by initialization script')
    w('=================================')
    w('')
    w('Most of :doc:`sph_input` applies to every run: the timestep controls, the')
    w('equation of state and the output intervals are read by the integrator and')
    w('mean the same thing whatever ``INAME`` is set to.')
    w('')
    w('The settings below are the exception.  Each is read only by the setup')
    w('routines listed, so setting one for any other script has no effect and')
    w('produces no warning.  This page is generated from the source, by looking at')
    w('which file actually mentions each variable.')
    w('')
    byscript = {}
    for v in members:
        codes, evo = usage.get(v.lower(), ([], False))
        if evo:
            continue
        for c in codes:
            byscript.setdefault(c, []).append(v)
    for c in sorted(byscript):
        title = '``%s``, %s' % (c, SCRIPT_TITLES.get(c, ''))
        w(title)
        w('-' * len(title))
        w('')
        w('.. list-table::')
        w('   :header-rows: 1')
        w('   :widths: 26 20 54')
        w('')
        w('   * - Variable'); w('     - Default'); w('     - Meaning')
        for v in sorted(byscript[c]):
            d, cm = defaults.get(v.lower(), ('--', ''))
            w('   * - ``%s``' % v)
            w('     - ``%s``' % d)
            w('     - %s' % (rst_escape(cm) if cm else '*Not yet documented.*'))
        w('')
    write_if_changed(out_path, '\n'.join(o) + '\n')
    return len(byscript)


def write_if_changed(path, text):
    """Write only when the content differs.

    These files live under source/, so rewriting one with identical content
    still updates its mtime, and a file watcher cannot tell that from a real
    edit.  With the generator running as a pre-build step that is a loop: build
    rewrites the file, the watcher sees a change, and it builds again.
    """
    try:
        if io.open(path, encoding='utf-8').read() == text:
            return False
    except IOError:
        pass
    io.open(path, 'w', encoding='utf-8').write(text)
    return True


def parse(path):
    lines = open(path, errors='ignore').read().split('\n')

    members, grabbing = [], False
    for l in lines:
        if re.match(r'\s+namelist/input/', l):
            grabbing = True
        elif grabbing and not re.match(r'\s{5}[$&]', l):
            grabbing = False
        if grabbing:
            l = re.sub(r'^\s+namelist/input/', '', re.sub(r'^\s{5}[$&]', '', l))
            members.append(l)
    members = [v.strip() for v in ','.join(members).split(',') if v.strip()]

    hi = next(i for i, l in enumerate(lines)
              if re.match(r"\s+open\(12,file='sph\.input'", l))
    lo = next(i for i in range(hi, 0, -1) if 'displacex=0d0' in lines[i])

    defaults = {}
    for l in lines[lo:hi]:
        if is_comment(l):
            continue
        m = re.match(r"\s+([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([^!]+?)\s*(?:!\s*(.*))?$", l)
        if m:
            defaults[m.group(1).lower()] = (m.group(2).strip(), (m.group(3) or '').strip())
    return members, defaults, lo + 1, hi

GROUPS = [
 ('Time and output', ['tf','dtout','nitpot','tjumpahead','tscanon','sepfinal','throwaway']),
 ('The particles',   ['n','nnopt','equalmass','starmass','starradius','hco','mco','hfloor']),
 ('Equation of state and physics',
                     ['neos','gam','nav','alpha','beta','ngr','nselfgravity','nkernel',
                      'ncooling','teq','reat']),
 ('Relaxation',      ['nrelax','trelax','treloff','tresplintmuoff','omega_spin']),
 ('Orbit of the encounter',
                     ['sep0','rp','impactparameter','e0','semimajoraxis','vinf2']),
 ('Compact object and black hole',
                     ['mbh','bbh_m1','bbh_m2','bbh_rp','bbh_semimajoraxis','bbh_vinf2',
                      'bbh_e0','bbh_trueanomaly','bbh_argperi','bbh_inclination',
                      'bbh_longitude']),
 ('Timestep control', ['nintvar','cn1','cn2','cn3','cn4','cn5','cn6','cn7']),
 ('Parallelism and GPUs',
                     ['ngravprocs','qthreads','ppn','computeexclusivemode','gflag']),
 ('Units',           ['runit','munit']),
 ('Input files',     ['startfile1','startfile2','startfile3','binaryfile','triplefile',
                      'bpbhfile','imagefile','advectedfile',
                      'eosfile','opacityfile','profilefile',
                      'stellarevolutioncodetype']),
]

def main():
    members, defaults, _lo, _hi = parse(INIT)
    usage = applies_to(os.path.dirname(INIT), members)
    placed = {v for _, vs in GROUPS for v in vs}
    leftover = [v for v in members if v.lower() not in placed]

    out = []
    w = out.append
    w('.. This page is generated by docs/generate_reference.py.  Do not edit by hand.')
    w('')
    w('sph.input')
    w('=========')
    w('')
    w('Every run is controlled by a Fortran namelist in ``sph.input``.  A minimal file')
    w('sets only what it needs, and anything omitted takes the default below.  The file')
    w('begins with ``&input`` and ends with ``&end``::')
    w('')
    w('    &input')
    w('     tf=200,')
    w('     n=20000,')
    w('    &end')
    w('')
    w('.. note::')
    w('')
    w('   The variables, defaults and descriptions on this page are read directly')
    w('   from the namelist declaration and the default-initialisation block in')
    w('   ``parallel_bleeding_edge/src/init.f``.  There are **%d** settings.'
      % len(members))
    w('')
    undoc = sum(1 for v in members if not defaults.get(v.lower(), ('', ''))[1])
    if undoc:
        w('.. warning::')
        w('')
        w('   %d of these have no description yet, because the corresponding line'
          % undoc)
        w('   in ``init.f`` carries no trailing comment.  Adding one there is all')
        w('   that is needed: it will appear here the next time this page is')
        w('   generated.')
        w('')
    w('')

    def scope(name):
        codes, evo = applies = usage.get(name.lower(), ([], False))
        if evo:
            return 'every run'
        if not codes:
            return '*nothing, see below*'
        return ' '.join('``%s``' % c for c in codes)

    def table(title, names):
        rows = [(n, defaults.get(n.lower(), ('--', ''))) for n in names
                if n.lower() in {m.lower() for m in members}]
        if not rows:
            return
        # A stable anchor per group, so other pages can link straight to it.
        slug = 'sph-input-' + re.sub(r'[^a-z0-9]+', '-', title.lower()).strip('-')
        w('.. _%s:' % slug); w('')
        w(title); w('-' * len(title)); w('')
        w('.. list-table::')
        w('   :header-rows: 1')
        w('   :widths: 19 16 20 45')
        w('')
        w('   * - Variable'); w('     - Default'); w('     - Applies to')
        w('     - Meaning')
        for n, (d, c) in rows:
            w('   * - ``%s``' % n)
            w('     - ``%s``' % d)
            w('     - %s' % scope(n))
            w('     - %s' % (rst_escape(c) if c else '*Not yet documented.*'))
        w('')

    for title, names in GROUPS:
        table(title, names)
    if leftover:
        table('Other settings', leftover)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    write_if_changed(OUT, '\n'.join(out) + '\n')
    n = write_by_script(os.path.join(os.path.dirname(OUT), 'settings_by_script.rst'),
                        members, defaults, usage)
    sys.stderr.write('wrote settings_by_script.rst (%d scripts)\n' % n)
    print('wrote %s  (%d settings, %d ungrouped)' % (OUT, len(members), len(leftover)))

if __name__ == '__main__':
    main()
