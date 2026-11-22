#!/usr/bin/env python3
"""sph_autotune.py: which ngravtree and ngpusph to write in sph.input, measured on this machine.

The automatic values (ngravtree=-1, ngpusph=-1) choose from the number of particles, the threads, the ranks and the GPUs, never from a time, so that a run and its
restart choose the same.  On hardware the rules were not fitted on, another pair of values can be faster.  This tool tries the pairs (ngravtree 0 or 2) x (ngpusph 0 or 7)
for a few tens of steps, on your input and on your machine, and says which two lines to write in sph.input.

    python3 sph_autotune.py [--dir RUNDIR] [--exe EXECUTABLE] [--launch "mpirun -np 1 --bind-to none"] [--threads T] [--steps 30] [--warmup 5]
                            [--reps 1] [--gravity 0,2 | keep] [--sph 0,7] [--work WORKDIR] [--keep] [--timeout SECONDS]

RUNDIR holds a run that is ready to start: sph.input, and sph.init (a new run) or restartrad.sph (a restart from a dump), with the input files.  Nothing in RUNDIR is changed.  Every candidate runs in a sub-folder of WORKDIR (default RUNDIR/autotune_<date>_<time>): small files are copied, large ones linked, and
sph.input is rewritten with ngravtree and ngpusph only.  The sub-folders are removed as soon as the candidate is done (--keep keeps them); autotune.txt and autotune.csv stay.
One candidate at a time; each is stopped (the whole process group: mpirun and the ranks) after warmup + steps steps.  A run with some steps much slower than the others (the machine
was busy with something else for a while) is measured again.

Standard library only; Python 3.6 or later."""
import argparse, csv, datetime, fnmatch, glob, os, re, select, shlex, shutil, signal, statistics, subprocess, sys, time

# files that the code writes in the folder it runs in (a link to one of the user's would be written through): copied, never linked
WRITTEN = ('restartrad.sph', 'm1m2rp.sph', 'parent.sph', 'startu.sph', 'polymc.col')
# what an earlier run left: not taken (the energy and log files are numbered from the first free number, so a copy would change the names the new run writes), unless sph.input
# or sph.init names the file (a start file can be an out*.sph)
LEFTOVER = ('out[0-9]*.sph', 'energy[0-9]*.sph', 'log[0-9]*.sph', 'col[0-9]*.sph', 'ecc*.sph', 'eat*.sph', 'jump*.sph', '*.old', 'core', 'core.[0-9]*')
COPY_LIMIT = 16 * 1024 * 1024        # bytes: a larger file is linked
DEFAULT_LAUNCH = 'mpirun -np 1 --bind-to none'    # one rank with every core free to the threads of OpenMP (the default of Open MPI for one rank is to bind it to one core)
EQUIVALENT = 0.03                    # two candidates within this fraction (or within the spread of the repetitions) are called equivalent
UNSURE = 10.0                        # percent: with one repetition a difference below this is not sure (a candidate has been seen to change by 10% between two runs)
UNSTEADY = 1.10                      # a run with the mean of its steps above this times their median is measured again (quiet runs gave 1.00 to 1.07; runs that shared the machine for a part of their steps 1.11 to 1.34)
TAIL = 6                            # lines of stderr shown for a candidate that died


class Usage(Exception):
    pass


class Interrupted(Exception):
    pass


# ------------------------------------------------------------------ sph.input ---

def split_comment(line):
    """(code, comment) of a namelist line: the comment starts at the first '!' that is not inside quotes."""
    q = None
    for i, ch in enumerate(line):
        if q:
            if ch == q:
                q = None
        elif ch in '\'"':
            q = ch
        elif ch == '!':
            return line[:i], line[i:]
    return line, ''


def find_end(code):
    """Index of the first '/' or '&end' of a line outside quotes (the end of the namelist group), or -1."""
    q = None
    for i, ch in enumerate(code):
        if q:
            if ch == q:
                q = None
        elif ch in '\'"':
            q = ch
        elif ch == '/':
            return i
        elif ch in '&$' and code[i + 1:i + 4].lower() == 'end':
            return i
    return -1


def key_regex(key):
    return re.compile(r'((?<![\w%%])%s\s*=\s*)([-+]?\d+)' % re.escape(key), re.I)


def nl_get(text, key):
    """The last integer written for key before the end of the namelist group, or None."""
    pat = key_regex(key)
    val = None
    for line in text.splitlines():
        code = split_comment(line)[0]
        end = find_end(code)
        for m in pat.finditer(code if end < 0 else code[:end]):
            val = int(m.group(2))
        if end >= 0:
            break
    return val


def nl_set(text, changes):
    """The text of a namelist with the keys of `changes` ({key: integer}) set: each value that is written is replaced (comments and the rest stay as they are), a key that is not
    written is added just before the end of the group."""
    nl = '\r\n' if '\r\n' in text else '\n'
    out = []
    done = set()
    ended = False
    for line in text.splitlines(True):
        body = line.rstrip('\r\n')
        eol = line[len(body):]
        if ended:
            out.append(line)
            continue
        code, comment = split_comment(body)
        end = find_end(code)
        head, tail = (code, '') if end < 0 else (code[:end], code[end:])
        for key, val in changes.items():
            head, n = key_regex(key).subn(lambda m, v=val: m.group(1) + str(v), head)
            if n:
                done.add(key)
        if end >= 0:
            ended = True
            todo = [k for k in changes if k not in done]
            if todo:
                if head.strip():
                    out.append(head.rstrip() + nl)
                    head = ''
                for k in todo:
                    out.append(' %s=%d,%s' % (k, changes[k], nl))
        out.append(head + tail + comment + eol)
    if not ended:
        raise Usage('sph.input has no end of the namelist group (a line with "&end" or "/")')
    return ''.join(out)


# --------------------------------------------------------------------- files ---

def prepare_dir(rundir, work, cdir, exe, input_text, changes, summary):
    """The folder of one candidate: the files of rundir copied (small) or linked (large), sph.input rewritten.  summary keeps the most bytes copied for one candidate and the names of the files linked."""
    os.makedirs(cdir)
    names = sorted(n for n in os.listdir(rundir) if not n.startswith('.'))
    try:
        init_text = open(os.path.join(rundir, 'sph.init'), errors='replace').read()
    except OSError:
        init_text = ''
    texts = [input_text, init_text]
    copied = 0
    for n in names:
        src = os.path.join(rundir, n)
        if n == 'sph.input' or n.startswith('autotune_') or os.path.realpath(src) in (os.path.realpath(exe), os.path.realpath(work)):
            continue
        if any(fnmatch.fnmatch(n, p) for p in LEFTOVER) and not any(re.search(r'(?<![\w.])%s(?![\w])' % re.escape(n), t) for t in texts):
            continue
        dst = os.path.join(cdir, n)
        if os.path.isdir(src):
            os.symlink(os.path.realpath(src), dst)
        elif os.path.isfile(src):
            size = os.path.getsize(src)
            if n in WRITTEN or size <= COPY_LIMIT:
                shutil.copyfile(src, dst)
                copied += size
            else:
                os.symlink(os.path.realpath(src), dst)
                summary['linked'].add(n)
    summary['copied'] = max(summary['copied'], copied)
    with open(os.path.join(cdir, 'sph.input'), 'w') as f:
        f.write(nl_set(input_text, changes))


def fingerprint(rundir, skip):
    """Name, size and modification time of every file of the folder (but not of the sub-folder `skip` nor of the autotune_* ones): what must be the same before and after."""
    fp = {}
    for root, dirs, files in os.walk(rundir):
        dirs[:] = [d for d in dirs if os.path.join(root, d) != skip and not d.startswith('autotune_')]
        for n in files:
            p = os.path.join(root, n)
            try:
                st = os.lstat(p)
            except OSError:
                continue
            fp[os.path.relpath(p, rundir)] = (st.st_size, st.st_mtime_ns)
    return fp


# ------------------------------------------------------------------ processes ---

def in_folder(cdir):
    """The pids of the processes whose working directory is cdir or below it (Linux, /proc; none elsewhere).  A candidate runs in a folder made for it, so these are its own: mpirun, the
    ranks, and anything they started that left the process group."""
    real = os.path.realpath(cdir)
    out = []
    try:
        names = os.listdir('/proc')
    except OSError:
        return out
    for p in names:
        if p.isdigit():
            try:
                cwd = os.readlink('/proc/%s/cwd' % p)
            except OSError:
                continue
            if cwd == real or cwd.startswith(real + os.sep):
                out.append(int(p))
    return out


def stop_tree(proc, cdir):
    """Stops a process started with start_new_session (the group: mpirun and the ranks) and every process that runs in cdir, and waits for them: SIGTERM, then SIGKILL."""
    if proc is None:
        return
    for sig, patience in ((signal.SIGTERM, 8.0), (signal.SIGKILL, 5.0)):
        try:
            os.killpg(proc.pid, sig)
        except OSError:
            pass
        for k in in_folder(cdir):
            try:
                os.kill(k, sig)
            except OSError:
                pass
        t0 = time.monotonic()
        while time.monotonic() - t0 < patience and (proc.poll() is None or in_folder(cdir)):
            time.sleep(0.05)
        if proc.poll() is not None and not in_folder(cdir):
            break
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass
    if proc.stdout:
        proc.stdout.close()


CURRENT = [None, None]       # the process and the folder being run: stopped by the handler of Ctrl-C, SIGTERM and SIGHUP
SIGNALLED = [False]


def on_signal(signum, frame):
    """The first signal interrupts the measurement (the cleaning that follows is not interrupted by a second one)."""
    if not SIGNALLED[0]:
        SIGNALLED[0] = True
        raise Interrupted('signal %d' % signum)


# ------------------------------------------------------------------- one run ---

AUTO_G = re.compile(r'gravity: ngravtree=-1 \(auto\):[^\n]*-> ngravtree=(\d)')
AUTO_S = re.compile(r'SPH on the GPU: ngpusph=-1 \(auto\):[^\n]*-> ngpusph=(\d+)')
SPH_PART = (('h, rho and neighbour lists', re.compile(r'^ngpusph: h, rho and neighbour lists (.*)$', re.M)),
            ('the pair sums (forces)', re.compile(r'^ngpusph: the pair sums \(forces\) (.*)$', re.M)))


def read_log(cdir):
    try:
        return open(os.path.join(cdir, 'log0.sph'), errors='replace').read()
    except OSError:
        return ''


def tail_lines(path, n=TAIL):
    try:
        return open(path, errors='replace').read().splitlines()[-n:]
    except OSError:
        return []


def reason_text(reason):
    return {'eof': 'the run ended by itself', 'timeout': 'time out'}.get(reason, str(reason))


class Runner(object):
    """Runs the candidates: one folder, one process group at a time."""

    def __init__(self, args, rundir, work, exe, input_text):
        self.args, self.rundir, self.work, self.exe, self.input_text = args, rundir, work, exe, input_text
        self.summary = dict(copied=0, linked=set())
        self.skip = max(args.warmup, 1)          # the first step has no step before it to take the difference with
        self.need = self.skip + args.steps       # the step lines to wait for
        self.count = 0
        self.info = {'ntot': None, 'openmp': None}

    def run(self, name, changes, probe=False):
        """One run with the namelist `changes`.  Returns a dict: status 'ok', 'died' or 'notmeasured'; ms (the times of the steps); note; log (log0.sph); stdout; stderr (its last lines);
        chose (probe only: the two values the automatic rule wrote, or None for the one it did not write); elapsed (seconds)."""
        args = self.args
        self.count += 1
        cdir = os.path.join(self.work, 'c%02d_%s' % (self.count, name))
        prepare_dir(self.rundir, self.work, cdir, self.exe, self.input_text, changes, self.summary)
        env = dict(os.environ)
        if args.threads:
            env['OMP_NUM_THREADS'] = str(args.threads)
        errf = open(os.path.join(cdir, 'stderr.log'), 'w')
        stamps, endit, lines = [], [], []   # arrival of the "mainit:" lines (the end of a step), of the " end it" lines of rank 0, and all the lines
        found = {'ngravtree': None, 'ngpusph': None}
        reason, proc, rc = None, None, None
        t0 = time.monotonic()
        try:
            try:
                proc = subprocess.Popen(args.launch_tokens + [self.exe], cwd=cdir, env=env, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=errf, start_new_session=True)
            except OSError as e:
                raise Usage('cannot start %s: %s' % (' '.join(args.launch_tokens + [self.exe]), e))
            CURRENT[0], CURRENT[1] = proc, cdir
            fd = proc.stdout.fileno()
            buf = b''
            deadline = t0 + args.timeout
            while reason is None:
                now = time.monotonic()
                if now >= deadline:
                    reason = 'timeout'
                    break
                if select.select([fd], [], [], min(0.2, deadline - now))[0]:
                    chunk = os.read(fd, 65536)
                    t = time.monotonic()
                    if not chunk:
                        reason = 'eof'
                        break
                    buf += chunk
                    while b'\n' in buf:
                        line, buf = buf.split(b'\n', 1)
                        text = line.decode('utf-8', 'replace')
                        if len(lines) < 20000:
                            lines.append(text)
                        if text.startswith('mainit:'):
                            stamps.append(t)
                        elif text.startswith(' end it'):
                            endit.append(t)
                if probe:
                    # stop when the lines of the choice are in log0.sph (or on the standard output, where the code writes them too)
                    text = read_log(cdir) + '\n' + '\n'.join(lines)
                    for key, pat in (('ngravtree', AUTO_G), ('ngpusph', AUTO_S)):
                        m = pat.search(text)
                        if m and found[key] is None:
                            found[key] = int(m.group(1))
                    if all(found[k] is not None for k in found if changes.get(k) == -1):
                        reason = 'probe'
                elif len(stamps) >= self.need or (not stamps and len(endit) >= self.need):
                    reason = 'enough'
            if reason == 'eof':
                try:
                    rc = proc.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    pass
        finally:
            stop_tree(proc, cdir)
            CURRENT[0] = CURRENT[1] = None
            errf.close()
        elapsed = time.monotonic() - t0
        err = tail_lines(os.path.join(cdir, 'stderr.log'))
        log = read_log(cdir)
        res = dict(status='ok', ms=[], note='', log=log, stdout='\n'.join(lines), stderr=err, rc=rc, elapsed=elapsed, chose=None)
        m = re.search(r'\bntot=\s*(\d+)', log)
        if m and self.info['ntot'] is None:
            self.info['ntot'] = int(m.group(1))
        m = re.search(r'OpenMP: P=\s*(\d+) T=\s*(\d+)', log)
        if m and self.info['openmp'] is None:
            self.info['openmp'] = (int(m.group(1)), int(m.group(2)))
        died = reason == 'eof' and rc not in (0, None)
        if probe:
            res['chose'] = (found['ngravtree'], found['ngpusph'])
            if reason != 'probe':
                res['status'] = 'died'
                res['note'] = self.why(reason, rc, err, len(stamps)) if died else 'the run ended without writing the choice of the automatic values (an older version of the code?)'
        else:
            marks = stamps if stamps else endit
            if died:
                res['status'], res['note'] = 'died', self.why(reason, rc, err, len(marks))
            elif len(marks) < self.skip + 3:
                res['status'] = 'notmeasured'
                res['note'] = '%s after %d steps in %.0f s, and at least %d are needed (%d of warm-up and 3 to time)' % (reason_text(reason), len(marks), elapsed, self.skip + 3, self.skip)
            else:
                res['ms'] = [(marks[k] - marks[k - 1]) * 1000.0 for k in range(self.skip, min(len(marks), self.need))]
                if sum(res['ms']) <= 0.0:
                    res['status'], res['ms'] = 'notmeasured', []
                    res['note'] = 'the lines of the %d steps arrived together (a busy machine, or an output held in a buffer), so the steps have no time' % len(marks)
                elif len(res['ms']) < args.steps:
                    res['note'] = '%s after %d steps: the %d steps after the warm-up are used' % (reason_text(reason), len(marks), len(res['ms']))
        if not args.keep and os.path.realpath(cdir).startswith(os.path.realpath(self.work) + os.sep):
            shutil.rmtree(cdir, ignore_errors=True)
        return res

    @staticmethod
    def why(reason, rc, err, steps):
        note = ('exit code %s after %d steps' % (rc, steps)) if reason == 'eof' else '%s after %d steps' % (reason_text(reason), steps)
        msg = [l.strip() for l in err if l.strip()]
        if msg:
            note += ': ' + msg[-1][:150]
        return note


def what_ran(text, s_asked):
    """What the log says was not done on the GPU although ngpusph asked for it (bit 1: the smoothing lengths, bit 2: the pair sums of the forces), as the line of log0.sph;
    '' when the log confirms it (or has no such line: a build that does not write them)."""
    bad = []
    for (what, pat), bit in zip(SPH_PART, (1, 2)):
        if s_asked & bit:
            found = pat.findall(text)
            if found and not found[-1].startswith('on the GPU'):
                bad.append('ngpusph: %s %s' % (what, found[-1].strip()))
    return '; '.join(bad)


# ----------------------------------------------------------------- statistics ---

class Cand(object):
    """A pair (ngravtree, ngpusph) and the runs made with it."""

    def __init__(self, g, s, kind='grid'):
        self.g, self.s, self.kind = g, s, kind
        self.runs = []
        self.status = 'notmeasured'
        self.note = ''
        self.ms = []
        self.rep_means = []

    @property
    def pair(self):
        return (self.g, self.s)

    def finish(self):
        """The status and the numbers from all the repetitions: a repetition that died or was changed by the code makes the candidate so."""
        bad = [r for r in self.runs if r['status'] != 'ok']
        order = {'died': 0, 'changed': 1, 'notmeasured': 2}
        if not self.runs:
            self.status = 'notrun'
        elif bad:
            bad.sort(key=lambda r: order[r['status']])
            self.status, self.note = bad[0]['status'], bad[0]['note']
        else:
            self.status = 'ok'
            self.note = '; '.join(sorted(set(r['note'] for r in self.runs if r['note'])))
        ok = self.status == 'ok'
        self.ms = [x for r in self.runs for x in r['ms']] if ok else []
        self.rep_means = [statistics.mean(r['ms']) for r in self.runs] if ok else []

    @property
    def steps(self):
        return len(self.ms) // max(len(self.runs), 1)

    @property
    def mean(self):
        return statistics.mean(self.ms) if self.ms else None

    @property
    def median(self):
        return statistics.median(self.ms) if self.ms else None

    @property
    def spread(self):
        """The distance between the largest and the smallest mean of the repetitions, relative to the smallest (None with one run)."""
        if len(self.rep_means) < 2:
            return None
        return (max(self.rep_means) - min(self.rep_means)) / min(self.rep_means)


def tolerance(a, b):
    """The fraction within which two candidates are equivalent: 3%, or the spread of the repetitions of either when that is larger."""
    return max([EQUIVALENT] + [c.spread for c in (a, b) if c.spread is not None])


def unsteady(ms):
    """True when some steps were so much slower than the others that the mean is above UNSTEADY times the median."""
    return statistics.mean(ms) > UNSTEADY * statistics.median(ms)


def measure(run_once, tell=lambda text: None):
    """The runs of one candidate that count, from `run_once()` (one run, as a dict of Runner.run): one, or two when the steps of the first were not steady.

    Something else that takes the machine for a few seconds makes some steps of one candidate slow, and that candidate loses for no reason of its own.  The steps of a run are the same
    steps every time, so a run is measured again when its steps were not steady: if the second run is steady the first is dropped (it was the machine), and if it is not steady
    either it is how this run goes (a step that writes a dump, lists that are built again) and both count."""
    r = run_once()
    if r['status'] != 'ok' or not unsteady(r['ms']):
        return [r]
    first = 'mean %.2f ms a step, median %.2f' % (statistics.mean(r['ms']), statistics.median(r['ms']))
    tell('the steps were not steady (%s): measured again' % first)
    again = run_once()
    if again['status'] != 'ok':
        return [r, again]
    if unsteady(again['ms']):
        again['note'] = '; '.join(x for x in (again['note'], 'the steps were not steady in two runs, so both count') if x)
        return [r, again]
    again['note'] = '; '.join(x for x in (again['note'], 'measured again: the steps of the first run were not steady (%s)' % first) if x)
    return [again]


def advise(cands, auto_pair, cur_pair=None):
    """(advised, fastest, equivalents) among the measured `cands`: the fastest; those within the tolerance of it; and from them, in this order of preference, the pair that sph.input gives
    now, the pair that the automatic rule chooses, the fastest."""
    ok = [c for c in cands if c.status == 'ok']
    if not ok:
        return None, None, []
    best = min(ok, key=lambda c: c.mean)
    eq = [c for c in ok if c is best or c.mean <= best.mean * (1 + tolerance(best, c))]
    rec = best
    for pair in (auto_pair, cur_pair):          # the later one wins
        for c in eq:
            if pair is not None and c.pair == pair:
                rec = c
    return rec, best, eq


def pct(x, ref):
    return 100.0 * (x - ref) / ref


# --------------------------------------------------------------------- output ---

STATUS_TEXT = {'ok': 'measured', 'died': 'died', 'changed': 'changed', 'notmeasured': 'not meas.', 'notrun': 'not run'}


def make_table(cands, best, marks):
    head = '%-24s %6s %9s %8s %8s %8s %8s %9s  %s' % ('candidate', 'steps', 'mean ms', 'median', 'min', 'max', 'vs best', 'status', 'note')
    rows = [head, '-' * len(head)]
    for c in cands:
        lab = 'ngravtree=%d ngpusph=%d' % c.pair
        spread = 'spread %.1f%% over %d repetitions' % (100 * c.spread, len(c.runs)) if c.spread is not None else ''
        note = '; '.join(x for x in (marks.get(c.pair, ''), spread, c.note) if x)
        if c.kind == 'extra':
            note = 'extra, not in the ranking' + ('; ' + note if note else '')
        if c.status == 'ok':
            vs = 'fastest' if c is best else '%+.1f%%' % pct(c.mean, best.mean)
            rows.append('%-24s %6d %9.2f %8.2f %8.2f %8.2f %8s %9s  %s' % (lab, c.steps, c.mean, c.median, min(c.ms), max(c.ms), vs, STATUS_TEXT['ok'], note))
        else:
            rows.append('%-24s %6s %9s %8s %8s %8s %8s %9s  %s' % (lab, '-', '-', '-', '-', '-', '-', STATUS_TEXT[c.status], note))
    return '\n'.join(rows)


CSV_HEADER = ['kind', 'ngravtree', 'ngpusph', 'chose_ngravtree', 'chose_ngpusph', 'status', 'steps', 'mean_ms', 'median_ms', 'min_ms', 'max_ms', 'percent_vs_fastest', 'spread_percent',
              'recommended', 'your_sph_input', 'note']


def csv_rows(cands, best, rec, rules, cur_pair):
    rows = []
    for c in cands:
        ok = c.status == 'ok'
        rows.append([c.kind, c.g, c.s, '', '', c.status, c.steps if ok else '', '%.3f' % c.mean if ok else '', '%.3f' % c.median if ok else '', '%.3f' % min(c.ms) if ok else '',
                     '%.3f' % max(c.ms) if ok else '', '%.2f' % pct(c.mean, best.mean) if ok and best else '', '%.2f' % (100 * c.spread) if c.spread is not None else '', 1 if c is rec else 0,
                     1 if c.kind == 'grid' and c.pair == cur_pair else 0, c.note])
    for kind, g, s, chosen, c, note in rules:
        ok = c is not None and c.status == 'ok'
        rows.append([kind, g, s, chosen[0] if chosen else '', chosen[1] if chosen else '', 'chose' if chosen else 'no_choice', '', '%.3f' % c.mean if ok else '', '%.3f' % c.median if ok else '',
                     '', '', '%.2f' % pct(c.mean, best.mean) if ok and best else '', '', 0, '', note])
    return rows


def advice_lines(cands, rec, best, eq, auto_pair, cur_pair, user_g, user_s, total_steps):
    """The advice as lines of text: the pair to write (or that nothing needs to change), what the choice changes, and what stays as it is.  cur_pair is the pair that sph.input gives now."""
    ranked = [c for c in cands if c.status == 'ok' and c.kind == 'grid']
    if rec is None:
        return ['No candidate was measured, so there is no advice.  The reasons are above.']
    cur = None
    for c in cands:
        if c.kind == 'grid' and cur_pair is not None and c.pair == cur_pair:
            cur = c
    has_auto = user_g == -1 or user_s == -1
    now = 'your sph.input as it is now (ngravtree=%d, ngpusph=%d%s)' % (user_g, user_s, ', which the code runs as ngravtree=%d, ngpusph=%d' % cur_pair if has_auto and cur_pair else '')
    out = ['Fastest: ngravtree=%d, ngpusph=%d (%.2f ms a step).' % (best.g, best.s, best.mean)]
    others = [c for c in eq if c is not best]
    reps = max(len(c.runs) for c in ranked)
    once = min(len(c.runs) for c in ranked) == 1          # some candidate has one run only: a difference of a few percent is not sure
    if others:
        out.append('Equivalent to it (within %.0f%%%s): %s.' % (100 * EQUIVALENT, ' or within the spread of the repetitions' if reps > 1 else '',
                                                                '; '.join('ngravtree=%d, ngpusph=%d (%+.1f%%)' % (c.g, c.s, pct(c.mean, best.mean)) for c in others)))
        if cur is not None and cur in eq:
            if cur is not best:
                out.append('The pair that your sph.input gives now is among them, so it is the one advised.')
        elif auto_pair is not None and any(c.pair == auto_pair for c in eq):
            out.append('The automatic rule (-1, -1) chooses ngravtree=%d, ngpusph=%d, which is among them, so that one is advised.' % auto_pair)
        else:
            out.append('Neither your sph.input as it is now nor the automatic rule (-1, -1) gives one of them, so the fastest is advised.')

    def block(c, title):
        return [title, '    ngravtree = %d' % c.g, '    ngpusph = %d' % c.s]

    def keep(title, ref):
        """Nothing to change: the pair of the file is the fastest or equivalent to it.  With a -1 in the file, the two lines that fix by hand what the code chooses now."""
        res = ['%s: %s gives %.2f ms a step, %s.' % (title, now, cur.mean, 'the fastest' if cur is ref[0] else 'within %.1f%% of %s' % (pct(cur.mean, ref[0].mean), ref[1]))]
        if has_auto:
            res += ['The file has -1, so the code chooses this pair itself; to fix it by hand (the choice can change on another machine or with another number of threads) write:',
                    '    ngravtree = %d' % cur.g, '    ngpusph = %d' % cur.s]
        return res

    tree = lambda g: g in (1, 2)
    main = rec
    refs = [(rec, eq, ranked)]                   # what each pair advised is measured against when the size of the differences is checked
    if user_g in (0, 1, 2) and tree(rec.g) != tree(user_g):
        same = [c for c in ranked if tree(c.g) == tree(user_g)]
        rec2, best2, eq2 = advise(same, auto_pair, cur_pair)
        out += ['', 'The fastest uses %s, and your sph.input has %s.  The tree is an approximation (treeacc), so the forces differ from the direct sum by the error of'
                % ('the tree' if tree(rec.g) else 'the direct sum', 'the tree' if tree(user_g) else 'the direct sum'),
                'the approximation, not only in the last digits.  The fastest with the gravity of your sph.input:']
        if rec2 is None:
            out.append('    (no candidate with that gravity was measured)')
        else:
            if rec2 is cur:
                out += keep('To keep the same gravity nothing needs to change', (best2, 'the fastest with that gravity'))
            else:
                out += block(rec2, 'To keep the same gravity write in sph.input:')
                if rec2.mean > rec.mean:
                    out.append('    (%.1f%% slower than the fastest overall)' % pct(rec2.mean, rec.mean))
            main = rec2
            refs.append((rec2, eq2, same))       # both pairs are advised here: the fastest overall and the fastest with the gravity of the file
        out += block(rec, 'To take the fastest and accept the other gravity solver write:')
    elif cur is not None and rec is cur:
        out += [''] + keep('Nothing to change', (best, 'the fastest'))
    else:
        out += [''] + block(rec, 'Write in sph.input:')
    # the pair that the file gives now: always said, with its time or why it has none
    if cur_pair is None or cur is None:
        out.append('Your %s: what the code chooses could not be read from log0.sph, so it was not measured.' % now[5:])
    elif cur.status != 'ok':
        out.append('Your %s: not measured (%s: %s).' % (now[5:], cur.status, cur.note))
    elif cur is not main:
        gain = 100.0 * (cur.mean - main.mean) / cur.mean
        if gain > 0:
            out.append('Your %s: %.2f ms a step; with the advised values %.2f, %.1f%% less.' % (now[5:], cur.mean, main.mean, gain))
        else:
            out.append('Your %s: %.2f ms a step, not slower than the advised values (%.2f).' % (now[5:], cur.mean, main.mean))
    # with one repetition a difference of a few percent is not sure: when the advice changes the values for a gain below UNSURE percent, or the nearest candidate that is not equivalent
    # is below UNSURE percent behind the one advised, the ranking can change from one run to the next
    if once:
        two = len(refs) == 2                     # two pairs are advised: the line says which difference it is about
        small = []
        if cur is not None and cur.status == 'ok' and cur is not main and cur.mean > main.mean:
            small.append((100.0 * (cur.mean - main.mean) / cur.mean, ' between your sph.input as it is now and the pair advised with the same gravity'))
        for ref, what in zip(refs, (' between the fastest overall and the next candidate', ' between the pair advised with the gravity of your sph.input and the next with that gravity')):
            slower = [c for c in ref[2] if c.status == 'ok' and c not in ref[1]]
            if slower:
                small.append((pct(min(c.mean for c in slower), ref[0].mean), what))
        small = [x for x in small if 0 < x[0] < UNSURE]
        if small:
            x, what = min(small)
            out += ['', 'The difference%s is %.1f%%: with one repetition a candidate can change by several percent from one run to the next; run again with --reps 3 to be sure.'
                    % (what if two else '', x)]
    # the SPH on the GPU changes the last digits of the forces
    if main is not cur and (cur_pair is None or (cur_pair[1] >= 3) != (main.s >= 3)):
        out += ['', 'ngpusph=0 and ngpusph=7 do not give the same last digits of the forces (see the page on the SPH on the GPUs): a run with the new value is not, bit for bit,',
                'the run with the old one.  Both are equally right.']
    out += ['',
            'The values written by hand stay as they are, so a run is reproducible.  sph_autotune changed nothing in your folder or in the code: write the lines yourself.',
            'The automatic values (-1) choose again from the particles, threads, ranks and GPUs, and may choose differently on another machine.',
            'This is the cost of the first %d steps of your run; measure again if the cost changes during the run (a collision), or if the machine, the threads or the particles change.' % total_steps]
    return out


# ----------------------------------------------------------------------- main ---

def parse_list(text, what):
    try:
        return [int(x) for x in text.replace(' ', '').split(',') if x != '']
    except ValueError:
        raise Usage('--%s takes a list of integers such as 0,2' % what)


def find_exe(rundir, given):
    if given:
        exe = os.path.abspath(given)
    else:
        found = sorted(set(glob.glob(os.path.join(rundir, '*_gpu_sph')) + glob.glob(os.path.join(rundir, '*_cpu_sph')) + glob.glob(os.path.join(rundir, 'run_sph'))))
        if len(found) != 1:
            raise Usage('no --exe given and %s in %s: say which executable to use' % (
                '%d candidates (%s)' % (len(found), ', '.join(os.path.basename(f) for f in found)) if found else 'no *_gpu_sph, *_cpu_sph or run_sph', rundir))
        exe = found[0]
    if not os.path.isfile(exe):
        raise Usage('the executable %s does not exist' % exe)
    if not os.access(exe, os.X_OK):
        raise Usage('the executable %s is not executable (chmod +x)' % exe)
    return exe


def say(text=''):
    sys.stdout.write(text + '\n')
    sys.stdout.flush()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', default='.', help='the folder of the run: sph.input, sph.init (a new run) or restartrad.sph (a restart), and the input files (default: the current folder)')
    ap.add_argument('--exe', help='the executable (default: the only *_gpu_sph, *_cpu_sph or run_sph of --dir)')
    ap.add_argument('--launch', default=DEFAULT_LAUNCH, help='the command that starts the executable, as in your own runs (default "%s"); the executable is added at the end; "" starts it directly'
                    % DEFAULT_LAUNCH)
    ap.add_argument('--threads', type=int, help='OpenMP threads per rank: sets OMP_NUM_THREADS for the runs (default: the environment as it is)')
    ap.add_argument('--steps', type=int, default=30, help='steps timed for each candidate, after the warm-up (default 30)')
    ap.add_argument('--warmup', type=int, default=5, help='steps before the timed ones (default 5)')
    ap.add_argument('--reps', type=int, default=1, help='repetitions of each candidate, in alternating order (default 1)')
    ap.add_argument('--gravity', help='ngravtree values to try, such as 0,2 (default 0,2 for a GPU build, 0,1 for a CPU build), or "keep": the value of sph.input stays')
    ap.add_argument('--sph', help='ngpusph values to try, such as 0,7 (default 0,7 for a GPU build, 0 for a CPU build)')
    ap.add_argument('--work', help='where the candidates run (default: RUNDIR/autotune_<date>_<time>); autotune.txt and autotune.csv are written there')
    ap.add_argument('--keep', action='store_true', help='keep the folders of the candidates')
    ap.add_argument('--timeout', type=float, default=600.0, help='seconds allowed to each run of a candidate (default 600)')
    args = ap.parse_args(argv)
    try:
        return tune(args)
    except Usage as e:
        sys.stderr.write('sph_autotune: %s\n' % e)
        return 2


def tune(args):
    t_start = time.monotonic()
    rundir = os.path.abspath(args.dir)
    if not os.path.isdir(rundir):
        raise Usage('the folder %s does not exist' % rundir)
    # a run that is ready to start has sph.input and either sph.init (a new run) or restartrad.sph (a restart from a dump: the code opens sph.init only when there is no restartrad.sph)
    if not os.path.isfile(os.path.join(rundir, 'sph.input')):
        raise Usage('sph.input not found in %s: the folder must hold a run that is ready to start' % rundir)
    if not os.path.isfile(os.path.join(rundir, 'sph.init')) and not os.path.isfile(os.path.join(rundir, 'restartrad.sph')):
        raise Usage('neither sph.init nor restartrad.sph found in %s: a run that is ready to start has sph.init (a new run) or restartrad.sph (a restart from a dump)' % rundir)
    if args.steps < 3 or args.warmup < 0 or args.reps < 1 or args.timeout <= 0:
        raise Usage('--steps must be at least 3, --warmup 0 or more, --reps 1 or more, --timeout above 0')
    exe = find_exe(rundir, args.exe)
    args.launch_tokens = shlex.split(args.launch) if args.launch.strip() else []
    if args.launch_tokens and not shutil.which(args.launch_tokens[0]):
        raise Usage('the launcher "%s" was not found: give --launch with the command you start your runs with, for example "mpirun -np 1 --bind-to none"' % args.launch_tokens[0])
    input_text = open(os.path.join(rundir, 'sph.input'), errors='replace').read()
    nl_set(input_text, {})                      # the namelist must have an end
    user_g = nl_get(input_text, 'ngravtree')
    user_g = 0 if user_g is None else user_g
    user_s = nl_get(input_text, 'ngpusph')
    user_s = -1 if user_s is None else user_s
    no_gravity = nl_get(input_text, 'ngr') == 0
    code_only = ''.join(split_comment(l)[0] + '\n' for l in input_text.splitlines())
    for m in re.finditer(r"['\"](\.\./[^'\"]*)['\"]", code_only):
        sys.stderr.write('sph_autotune: warning: sph.input names %s, above the folder of the run: the candidates run one folder deeper and will not find it\n' % m.group(1))
    work = os.path.abspath(args.work) if args.work else os.path.join(rundir, 'autotune_' + datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
    if os.path.isdir(work) and os.listdir(work):
        raise Usage('the work folder %s exists and is not empty' % work)
    try:
        os.makedirs(work, exist_ok=True)
    except OSError as e:
        raise Usage('cannot make the work folder %s: %s' % (work, e))
    before = fingerprint(rundir, work)
    runner = Runner(args, rundir, work, exe, input_text)
    SIGNALLED[0] = False
    old_handlers = [(sig, signal.signal(sig, on_signal)) for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)]
    interrupted = False
    say('sph_autotune: folder %s' % rundir)
    say('  executable %s' % exe)
    say('  launch "%s <executable>"%s' % (args.launch, ', OMP_NUM_THREADS=%d' % args.threads if args.threads else ''))
    say('  %d steps after %d of warm-up, %d repetition(s), work folder %s' % (args.steps, runner.skip, args.reps, work))
    say()
    cands, probes = [], {}
    cpu_build = exe.endswith('_cpu_sph')
    try:
        # what the automatic rule chooses: both values -1, the default case (the gravity of sph.input, ngpusph=-1), and, when sph.input has ngravtree=-1 and an ngpusph written, the file as it is
        for name, g, sv in (('auto', -1, -1), ('default', user_g, -1), ('yours', -1, user_s)):
            if name == 'default' and user_g == -1:
                probes[name] = probes['auto']
                continue
            if name == 'yours' and not (user_g == -1 and user_s >= 0):
                continue
            say('[%s] ngravtree=%d ngpusph=%d: only to read what the code chooses ...' % (name, g, sv))
            r = runner.run(name, {'ngravtree': g, 'ngpusph': sv}, probe=True)
            probes[name] = r
            say('      %s' % ('chose ' + ' '.join('%s=%d' % (k, v) for k, v in zip(('ngravtree', 'ngpusph'), r['chose']) if v is not None) if r['status'] == 'ok' else r['note']))
        m = AUTO_S.search(probes['auto']['log'] + '\n' + probes['auto']['stdout'])
        cpu_build = cpu_build or (m is not None and 'no GPU code in this build or no GPU' in m.group(0))
        if args.gravity is not None and args.gravity.strip().lower() == 'keep':
            glist = [user_g]
        else:
            glist = parse_list(args.gravity, 'gravity') if args.gravity else ([0, 1] if cpu_build else [0, 2])
        slist = parse_list(args.sph, 'sph') if args.sph else ([0] if cpu_build else [0, 7])
        if no_gravity and len(glist) > 1:
            glist = [user_g]
            say('ngr=0 in sph.input, no gravity: ngravtree is not tried (it stays %d)' % user_g)
        cands = [Cand(g, s) for g in glist for s in slist]
        # the pair that sph.input gives now is always measured, and ranked like the others
        cur_pair = current_pair(probes, user_g, user_s)
        if cur_pair is not None and cur_pair[0] >= 0 and cur_pair[1] >= 0 and cur_pair not in [c.pair for c in cands]:
            cands.append(Cand(cur_pair[0], cur_pair[1]))
        # the pair that a rule chooses, if it is not in the list, is measured as well: the table says how far it is from the fastest (it does not enter the ranking)
        for name in ('auto', 'default'):
            p = chosen_pair(probes[name], user_g, user_s)
            if p is not None and p[0] >= 0 and p[1] >= 0 and p not in [c.pair for c in cands]:
                cands.append(Cand(p[0], p[1], 'extra'))
        say('%s build; candidates: %s' % ('CPU' if cpu_build else 'GPU', ', '.join('(%d,%d)' % c.pair for c in cands)))
        say()
        total, k = len(cands) * args.reps, 0
        for rep in range(args.reps):
            for c in (cands if rep % 2 == 0 else cands[::-1]):
                k += 1
                say('[%d/%d] ngravtree=%d ngpusph=%d%s ...' % (k, total, c.g, c.s, ' (repetition %d)' % (rep + 1) if args.reps > 1 else ''))

                def run_once(c=c):
                    r = runner.run('g%d_s%d' % c.pair, {'ngravtree': c.g, 'ngpusph': c.s})
                    if r['status'] == 'ok' and c.s >= 1:
                        bad = what_ran(r['log'] + '\n' + r['stdout'], c.s)
                        if bad:
                            r['status'], r['note'] = 'changed', 'log0.sph says "%s"' % bad
                    return r

                runs = measure(run_once, lambda text: say('      ' + text))
                c.runs += runs
                r = runs[-1]
                if r['status'] == 'ok':
                    say('      %.2f ms a step (median %.2f, %d steps, %.0f s)%s' % (statistics.mean(r['ms']), statistics.median(r['ms']), len(r['ms']), r['elapsed'], '; ' + r['note'] if r['note'] else ''))
                else:
                    say('      %s: %s' % (r['status'], r['note']))
    except Usage:
        try:
            os.rmdir(work)
        except OSError:
            pass
        raise
    except (KeyboardInterrupt, Interrupted):
        interrupted = True
        say()
        say('sph_autotune: interrupted; the run in progress is stopped and the folders of the candidates removed')
    finally:
        if CURRENT[0] is not None:
            stop_tree(CURRENT[0], CURRENT[1])
        for sig, h in old_handlers:
            signal.signal(sig, h)
        if not args.keep:
            for d in glob.glob(os.path.join(work, 'c[0-9][0-9]_*')):
                if os.path.realpath(d).startswith(os.path.realpath(work) + os.sep):
                    shutil.rmtree(d, ignore_errors=True)
    for c in cands:
        c.finish()
    after = fingerprint(rundir, work)
    changed = sorted(p for p in set(before) | set(after) if before.get(p) != after.get(p))
    text, rec, best, rules = report(args, runner, rundir, cands, probes, user_g, user_s, interrupted, changed, len(before), time.monotonic() - t_start)
    say()
    say(text)
    try:
        with open(os.path.join(work, 'autotune.txt'), 'w') as f:
            f.write(text + '\n')
        with open(os.path.join(work, 'autotune.csv'), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(CSV_HEADER)
            w.writerows(csv_rows(cands, best, rec, rules, current_pair(probes, user_g, user_s)))
        say()
        say('Written: %s and %s' % (os.path.join(work, 'autotune.txt'), os.path.join(work, 'autotune.csv')))
    except OSError as e:
        sys.stderr.write('sph_autotune: could not write the table: %s\n' % e)
    if interrupted:
        return 130
    return 0 if rec is not None and not changed else 1


def chosen_pair(probe, user_g, user_s):
    """The (ngravtree, ngpusph) that a run with the automatic values would use: the values the code wrote, and the user's for the one that was not -1; None if the log says nothing."""
    if probe is None or probe['status'] != 'ok':
        return None
    g, s = probe['chose']
    return (g if g is not None else user_g, s if s is not None else user_s)


def current_pair(probes, user_g, user_s):
    """The (ngravtree, ngpusph), with no -1, that the sph.input of the user gives now: both values as written, or the -1 completed with what the code chose in the probe (the gravity of the
    file with ngpusph=-1; both -1; ngravtree=-1 with the ngpusph of the file).  None if the log said nothing."""
    if user_g >= 0 and user_s >= 0:
        return (user_g, user_s)
    return chosen_pair(probes.get('default' if user_g >= 0 else ('auto' if user_s < 0 else 'yours')), user_g, user_s)


def report(args, runner, rundir, cands, probes, user_g, user_s, interrupted, changed, nfiles, seconds):
    """The text of the result and what the CSV needs: (text, advised candidate, fastest, rows of the rules)."""
    ranked = [c for c in cands if c.kind == 'grid']
    best = None
    ok = [c for c in ranked if c.status == 'ok']
    if ok:
        best = min(ok, key=lambda c: c.mean)
    by_pair = dict((c.pair, c) for c in cands)
    auto_pair = chosen_pair(probes.get('auto'), user_g, user_s)
    def_pair = chosen_pair(probes.get('default'), user_g, user_s)
    cur_pair = current_pair(probes, user_g, user_s)
    marks = {}
    for p, label in ((cur_pair, 'your sph.input'), (auto_pair, 'chosen by the automatic values'), (def_pair if user_g != -1 else None, 'chosen by ngpusph=-1 with your gravity')):
        if p is not None and p in by_pair:
            marks[p] = (marks[p] + ', ' if p in marks else '') + label
    marks = dict((p, '<- ' + v) for p, v in marks.items())
    rec, fastest, eq = advise(ranked, auto_pair, cur_pair) if not interrupted else (None, None, [])
    info = runner.info
    lines = ['Result of sph_autotune%s for %s%s%s' % (' (INTERRUPTED: only what was measured)' if interrupted else '', rundir, ', %d particles' % info['ntot'] if info['ntot'] else '',
                                                     ', %d rank(s) with %d thread(s) (from log0.sph)' % info['openmp'] if info['openmp'] else ''), '']
    lines.append(make_table(cands, best, marks) if cands else '(no candidate was run)')
    lines.append('')
    rules = []
    for name, label, g, sv in (('auto', 'The automatic values (ngravtree=-1, ngpusph=-1)', -1, -1),
                               ('default', 'The default case (ngravtree=%d as in your sph.input, ngpusph=-1)' % user_g, user_g, -1),
                               ('yours', 'Your sph.input (ngravtree=-1, ngpusph=%d)' % user_s, -1, user_s)):
        if probes.get(name) is None or (name == 'default' and user_g == -1):
            continue
        p = {'auto': auto_pair, 'default': def_pair, 'yours': cur_pair}[name]
        if p is None:
            note = probes[name]['note'] or 'log0.sph has no line with the choice (an older version of the code?)'
            lines.append('%s: %s' % (label, note))
            rules.append((name, g, sv, None, None, note))
            continue
        c = by_pair.get(p)
        if c is not None and c.status == 'ok' and best is not None:
            if c is best:
                how = 'the fastest measured'
            elif c.mean < best.mean:
                how = '%.2f ms, %.1f%% faster than the fastest of the ranking' % (c.mean, -pct(c.mean, best.mean))
            else:
                how = '%.2f ms against %.2f for the fastest, %.1f%% slower' % (c.mean, best.mean, pct(c.mean, best.mean))
            lines.append('%s: the code chooses ngravtree=%d ngpusph=%d: %s.' % (label, p[0], p[1], how))
        else:
            lines.append('%s: the code chooses ngravtree=%d ngpusph=%d, which was not measured%s.' % (label, p[0], p[1], ' (%s)' % c.note if c is not None and c.note else ''))
        rules.append((name, g, sv, p, c, ''))
    lines.append('')
    for c in cands:
        if c.status not in ('ok', 'notrun'):
            lines.append('ngravtree=%d ngpusph=%d: %s (%s).' % (c.g, c.s, {'died': 'the run died', 'changed': 'not run as asked', 'notmeasured': 'not measured'}[c.status], c.note))
            for r in c.runs:
                if r['status'] == 'died':
                    lines += ['      stderr: %s' % l.rstrip() for l in r['stderr'][-TAIL:]]
                    break
            lines.append('      It is left out of the ranking.')
    if any(c.status not in ('ok', 'notrun') for c in cands):
        lines.append('')
    if not interrupted:
        lines += advice_lines(cands, rec, best, eq, auto_pair, cur_pair, user_g, user_s, runner.skip + args.steps)
    lines.append('')
    if changed:
        lines.append('WARNING: files of your folder changed during the measurement: %s' % ', '.join(changed[:10]))
    else:
        lines.append('Your folder is unchanged (name, size and modification time of %d files, before and after).' % nfiles)
    lines.append('Each candidate had %.1f MB copied and %d large files linked; the folders of the candidates %s.  Time of sph_autotune: %.0f s.' % (
        runner.summary['copied'] / 1048576.0, len(runner.summary['linked']), 'are kept' if args.keep else 'were removed', seconds))
    return '\n'.join(lines), rec, best, rules


if __name__ == '__main__':
    sys.exit(main())
