#!/usr/bin/env python3
"""
run_all.py -- the weekly gap-driven benchmark job.

    make bench-gap                       # everything, three systems
    python3 benchmarks/run_all.py --list
    python3 benchmarks/run_all.py --only 01,31 --timeout 60
    python3 benchmarks/run_all.py --system mathilda,python

Runs every `benchmarks/NN-slug/` experiment in Mathilda, in Python
(numpy/scipy/sympy), and -- when `wolframscript` is on the box -- in
Mathematica.  Joins the cases by label, classifies each one, and writes:

    REPORT.md                  ranked SLOWER cases + the coverage percentage
    ABSENT.md                  the absence register: what we do not have
    results/<YYYY-MM-DD>.json  every raw case plus the host/version block

WHY A JOB AND NOT ANOTHER WRITE-UP.  `docs/experiments/` holds twenty
experiments as narrative, and `docs/experiments/run.sh` only prints the three
systems' raw text side by side -- it does not parse, join, rank, or classify.
`comparisons/hpc_bench.py` does aggregate, but its kernels are inline in a
1791-line Python file, so there is nothing on disk to edit or run standalone.
This runs from kept `.m`/`.py` pairs and produces a ranked report, so it can be
re-run at the end of every week and the diff read as the week's dev work.

THE FIVE CLASSES, and why absence is never a ratio
--------------------------------------------------
    SLOWER      every system answered; Mathilda is >= SLOWER_AT x the best other
    AHEAD       Mathilda is fastest, or within SLOWER_AT x of the best
    ABSENT      the head does not exist in Mathilda -- NO ratio, ever
    INCOMPLETE  the head exists but returned unevaluated, errored, or timed out
    CHECK-FAIL  the systems disagree on the check value -- timing is DISCARDED

`SLOWER` is kernel work and has a ratio.  `ABSENT` is feature work and must not:
pooling "we do not have DSolve" into a speed average turns a missing algorithm
into a fake 40x.  A `CHECK-FAIL` case prints no timing at all -- an unverified
timing is worse than a missing one, because the columns may be timing different
programs.  Research found exactly that: all 15 WolframMark tests spell their
data `RandomReal[{}, dims]`, which Mathematica accepts and Mathilda returns
UNEVALUATED, so eight of them "ran" in ~0.1 ms while computing nothing.

WHAT IS NOT TRUSTED
-------------------
The exit code.  Mathilda exits 0 after `Power::infy` and after an arity error,
with diagnostics on stderr.  Classification comes from the parsed cases plus
captured stderr, never from `$?`.
"""

import argparse
import datetime
import json
import os
import platform
import re
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BENCH = os.path.join(ROOT, "benchmarks")

# Tool locations.  The defaults in the older harnesses (a hardcoded
# /usr/local/bin/python3.11 and an /Applications/... wolframscript) are both
# wrong on at least one machine this has run on, which is why every one is an
# env override with an autodetected fallback.
MATHILDA = os.environ.get("MATHILDA_BIN", os.path.join(ROOT, "Mathilda"))
PYTHON = os.environ.get("HPC_PYTHON", sys.executable or "python3")


def find_wolfram():
    env = os.environ.get("WOLFRAMSCRIPT")
    if env:
        return env if os.path.exists(env) else None
    onpath = shutil.which("wolframscript")
    if onpath:
        return onpath
    mac = "/Applications/Mathematica.app/Contents/MacOS/wolframscript"
    return mac if os.path.exists(mac) else None


WOLFRAM = find_wolfram()

# A case is a finding once it is at least this many times the best other system.
# 1.5 rather than 1.0 because run-to-run spread on large-working-set cases is
# documented at +-20% in docs/experiments/README.md; 1.5 sits outside that.
SLOWER_AT = 1.5

# Check-value agreement.  Floats compare on relative error; anything else is
# compared as normalised text.
CHECK_RTOL = 1e-6

SYSTEMS = ("mathilda", "wolfram", "python")
PRETTY = {"mathilda": "Mathilda", "wolfram": "Mathematica", "python": "Python"}


# --------------------------------------------------------------------------
# Discovery
# --------------------------------------------------------------------------

EXP_RE = re.compile(r"^(\d{2})-(.+)$")


def discover(only=None):
    """Every benchmarks/NN-slug/ folder holding at least one .m or .py."""
    out = []
    for name in sorted(os.listdir(BENCH)):
        m = EXP_RE.match(name)
        if not m:
            continue
        d = os.path.join(BENCH, name)
        if not os.path.isdir(d):
            continue
        num, slug = m.group(1), m.group(2)
        if only and num not in only and slug not in only and name not in only:
            continue
        mfile = next((f for f in sorted(os.listdir(d)) if f.endswith(".m")), None)
        pyfile = next((f for f in sorted(os.listdir(d)) if f.endswith(".py")), None)
        if not mfile and not pyfile:
            continue
        out.append({"num": num, "slug": slug, "name": name, "dir": d,
                    "m": mfile, "py": pyfile,
                    "group": group_of(int(num))})
    return out


def group_of(num):
    """Experiment number -> area. 01-10 / 11-20 / 21-30 by construction.

    E (53-62) is advanced numerical analysis against scipy/numpy (one mpmath
    baseline) — decompositions, eigenproblems, multidim/oscillatory quadrature,
    stiff ODE/PDE, nonlinear systems, high-degree roots, DCT/DST, regularized
    least squares, arbitrary precision. It sits one level below group B on the
    same subsystems and is run under Python 3.11 (see benchmarks/requirements.txt).


    C is exactly the eight open roadmap items from docs/experiments/README.md —
    the array substrate: packing, dispatch, per-call overhead. D is the two
    subsystems (src/graph/, src/strings/) that nothing in the corpus had ever
    timed; they are neither roadmap items nor numeric, so folding them into C
    would have mislabelled both.

    There was a different group D holding Wolfram's WolframMark suite. It was dropped:
    WolframMark is a HARDWARE benchmark — it holds Mathematica constant and
    varies the machine, scoring against a reference system — so its sizes encode
    Mathematica's performance profile rather than the workloads' intrinsic cost,
    and it measures nothing about a CAS's quality or coverage. Its one durable
    finding (RandomReal[{}, dims] returns unevaluated in Mathilda where
    Mathematica reads it as the default {0,1} range) is recorded in the changelog;
    a guard for it belongs in tests/, not in a timing job."""
    if 1 <= num <= 10:
        return "A symbolic (sympy)"
    if 11 <= num <= 20:
        return "B numeric libraries (scipy)"
    if 21 <= num <= 28:
        return "C array substrate (numpy)"
    if 53 <= num <= 62:
        return "E numerical analysis (scipy)"
    return "D uncovered subsystems"


# --------------------------------------------------------------------------
# Running one file, in one system
# --------------------------------------------------------------------------

LINE_RE = re.compile(r"^(BENCH|CHECK|REQUIRE|SKIP)\t([^\t]*)\t(.*)$")


def fmt_dur(s):
    """Human durations: 42s, 3m07s, 1h02m."""
    s = max(0, int(round(s)))
    if s < 60:
        return "%ds" % s
    m, sec = divmod(s, 60)
    if m < 60:
        return "%dm%02ds" % (m, sec)
    h, m = divmod(m, 60)
    return "%dh%02dm" % (h, m)


BAR_W = 26

# How often the poll loop redraws. 0.25 s is smooth to the eye and costs nothing
# next to subprocesses that run for seconds to minutes.
TICK_S = 0.25

# When output is NOT a tty (piped to tee, redirected to a log, captured by CI)
# carriage returns are useless -- they would collapse into one unreadable line.
# Instead emit a fresh heartbeat line this often, so a long experiment still
# shows movement in the log rather than appearing hung for four minutes.
HEARTBEAT_S = 15.0


def bar(frac):
    """A 26-cell block bar. frac is clamped to [0, 1]."""
    frac = min(1.0, max(0.0, frac))
    full = int(frac * BAR_W)
    return "\u2588" * full + "\u2591" * (BAR_W - full)


class Progress:
    """Draws the one-line bar, and keeps drawing it while a subprocess runs.

    Two output modes, chosen by whether stderr is a terminal:

      tty      redraw in place with \r every TICK_S -- a smoothly moving bar.
      not tty  a fresh line every HEARTBEAT_S. Carriage returns in a log file
               collapse into a single unreadable line, but a log that shows
               nothing for the four minutes of a timeout experiment looks hung,
               so it still has to move -- just more slowly, and by appending.

    The bar is weighted by TIME, not by experiment count, whenever a previous
    run is available to weight with. Counting experiments would show 87% while
    the two four-minute timeout experiments were still to come.
    """

    def __init__(self, total, exps, prev, verbose):
        self.total, self.exps, self.prev = total, exps, prev
        self.verbose = verbose
        self.tty = sys.stderr.isatty()
        self.t_start = None
        self.i, self.slug = 0, ""
        self.done_names = []
        self.settled = 0.0        # elapsed as of the last COMPLETED experiment
        self.t_exp = None         # when the in-flight experiment began
        self.cap = None           # most one unfinished experiment can still owe
        self.sys_label = ""
        self.last_emit = 0.0
        self.dirty = False

    def start(self, t_start):
        self.t_start = t_start

    def begin(self, i, slug):
        """A new experiment: always visible, in either mode."""
        import time as _t
        self.i, self.slug = i, slug
        self.sys_label = ""
        self.t_exp = _t.time()
        self.draw(boundary=True)

    def system(self, name):
        """Which of the three systems is running. Redraws on a tty; in a log it
        rides the heartbeat, since three extra lines per experiment would treble
        the noise for no information a heartbeat does not already carry."""
        self.sys_label = name
        self.draw()

    def finish_experiment(self, name):
        import time as _t
        self.done_names.append(name)
        self.settled = _t.time() - self.t_start   # ratio basis: completed only
        self.t_exp = None

    def _frac_and_eta(self):
        import time as _t
        elapsed = _t.time() - self.t_start
        remaining = [x["name"] for x in self.exps[self.i - 1:]] if self.i else \
            [x["name"] for x in self.exps]
        in_flight = (_t.time() - self.t_exp) if self.t_exp else None
        rem = eta(self.done_names, remaining, elapsed, self.prev,
                  settled=self.settled, in_flight=in_flight, cap=self.cap)
        if self.prev:
            tot = sum(self.prev.get(e["name"], 0) for e in self.exps) or 1
            done = sum(self.prev.get(n, 0) for n in self.done_names)
            # Creep through the in-flight experiment too, so the bar keeps
            # moving during a single long one instead of stepping at the end.
            cur = self.prev.get(self.exps[self.i - 1]["name"], 0) if self.i else 0
            spent_here = elapsed - sum(self.prev.get(n, 0) for n in self.done_names) \
                * (elapsed / done if done else 1)
            partial = min(cur, max(0.0, spent_here)) if cur else 0.0
            frac = (done + partial) / tot
        else:
            frac = (self.i - 1) / float(self.total) if self.total else 0.0
        return min(0.999, max(0.0, frac)), elapsed, rem

    def draw(self, boundary=False, final=False):
        import time as _t
        if self.t_start is None:
            return
        now = _t.time()
        if final:
            frac, elapsed, rem = 1.0, now - self.t_start, 0
        else:
            frac, elapsed, rem = self._frac_and_eta()

        line = "  %s %3d%%  %2d/%d  %-7s  %-9s %s%s" % (
            bar(frac), round(100 * frac), self.i, self.total, fmt_dur(elapsed),
            ("ETA " + fmt_dur(rem)) if rem is not None else "ETA --",
            (self.slug if not final else "done")[:24],
            ("  [%s]" % self.sys_label) if self.sys_label and not final else "")

        if self.tty and not final:
            sys.stderr.write("\r" + line + " " * 8)
            sys.stderr.flush()
        else:
            # An experiment boundary always prints; ticks in between are
            # rate-limited so a log shows movement without filling with it.
            if not (boundary or final) and (now - self.last_emit) < HEARTBEAT_S:
                return
            self.last_emit = now
            print(line, file=sys.stderr, flush=True)

    def note(self, text):
        """A line that must survive the bar (a timeout, a per-system count)."""
        if self.tty:
            sys.stderr.write("\r" + " " * 78 + "\r")
        print(text, file=sys.stderr, flush=True)
        self.last_emit = 0.0

    def close(self):
        self.draw(final=True)
        if self.tty:
            sys.stderr.write("\n")
            sys.stderr.flush()


def load_prev_durations():
    """Per-experiment wall times from the most recent previous run, if any.

    This is what makes the ETA worth printing.  Experiment costs here span three
    orders of magnitude (a 0.3 ms Expand case next to a 23 s Fit case), so a naive
    "elapsed / completed * remaining" extrapolation is wildly wrong early on --
    it will happily predict 4 minutes while experiment 29 alone is about to take
    six.  Replaying the last run's actual per-experiment times, rescaled by how
    this run is tracking against them so far, is accurate within a few percent
    from the second experiment onward.
    """
    outdir = os.path.join(BENCH, "results")
    if not os.path.isdir(outdir):
        return {}
    cands = [os.path.join(outdir, n) for n in os.listdir(outdir)
             if n.endswith(".json")]
    for path in sorted(cands, key=lambda q: -os.path.getmtime(q)):
        try:
            with open(path) as f:
                d = json.load(f)
            if d.get("durations") and not d.get("partial"):
                return d["durations"]
        except Exception:
            continue
    return {}


def eta(done_names, remaining_names, elapsed, prev, settled=None,
        in_flight=None, cap=None):
    """Seconds remaining. Replays `prev`, rescaled by how this run is tracking.

    `settled` is the elapsed time as of the last COMPLETED experiment, and it is
    what the tracking ratio is computed from. Using total `elapsed` instead lets
    time spent inside the in-flight experiment inflate the ratio, which is
    catastrophic during a timeout: stuck four minutes into `series-limit`, the
    ratio tripled and the ETA climbed 25m -> 1h53m while the true remainder was
    unchanged. A stall must shrink the current experiment's own estimate toward
    zero, never re-scale the twenty experiments behind it.

    `in_flight` is time already spent in the current experiment, subtracted from
    that experiment's own predicted cost (floored at zero). `cap` bounds what a
    single unfinished experiment can still owe — a hung one cannot exceed its
    timeout budget, so the ETA stops growing rather than running away.
    """
    if settled is None:
        settled = elapsed
    known = [n for n in remaining_names if n in prev]
    unknown = [n for n in remaining_names if n not in prev]

    if done_names:
        prev_done = sum(prev.get(n, 0) for n in done_names)
        ratio = (settled / prev_done) if prev_done > 0 else 1.0
        # A single slow experiment should not swing the whole projection.
        ratio = min(3.0, max(0.33, ratio))
    else:
        ratio = 1.0

    if not known and not unknown:
        return 0.0

    total = 0.0
    for idx, n in enumerate(remaining_names):
        pred = prev.get(n)
        if pred is None:
            pred = (settled / len(done_names)) if done_names else 0.0
        else:
            pred *= ratio
        if idx == 0 and in_flight is not None:
            pred = max(0.0, pred - in_flight)
            if cap is not None:
                pred = min(pred, max(0.0, cap - in_flight))
        total += pred
    return total


def run_one(system, exp, timeout, tick=None):
    """Launch one file and parse its tagged lines.

    cwd is the experiment folder, which is what makes `Get["../harness.m"]`
    resolve -- Mathilda's Get is relative to cwd and $InputFileName does not
    exist.  Partial output survives a timeout kill because repl.c line-buffers
    stdout when a script's output is redirected.
    """
    if system == "mathilda":
        if not exp["m"] or not os.path.exists(MATHILDA):
            return None
        cmd = [MATHILDA, "-file", exp["m"]]
    elif system == "wolfram":
        if not exp["m"] or not WOLFRAM:
            return None
        cmd = [WOLFRAM, "-file", exp["m"]]
    else:
        if not exp["py"]:
            return None
        cmd = [PYTHON, exp["py"]]

    # Poll rather than block, so the caller can animate while this runs.
    #
    # subprocess.run() blocks until exit, which froze the progress bar for the
    # whole of a 240 s experiment -- the display only moved at experiment
    # boundaries, which is exactly when you least need it. Polling lets `tick`
    # redraw a few times a second.
    #
    # Output goes to TEMPORARY FILES rather than PIPE. With PIPE, a child that
    # fills the 64 KiB pipe buffer blocks forever unless someone drains it, and
    # draining it concurrently means threads. Files cannot deadlock, and these
    # children already write line-buffered output we only read at the end.
    import tempfile
    import time as _time

    timed_out = False
    fo = tempfile.TemporaryFile(mode="w+", encoding="utf-8", errors="replace")
    fe = tempfile.TemporaryFile(mode="w+", encoding="utf-8", errors="replace")
    try:
        proc = subprocess.Popen(cmd, cwd=exp["dir"], stdout=fo, stderr=fe,
                                text=True)
    except OSError as e:
        fo.close(); fe.close()
        return {"error": str(e), "bench": {}, "check": {}, "require": {},
                "skip": {}, "stderr": "", "timed_out": False, "rc": None}

    t0 = _time.time()
    while True:
        rc = proc.poll()
        if rc is not None:
            break
        if _time.time() - t0 > timeout:
            proc.kill()
            proc.wait()
            rc, timed_out = None, True
            break
        if tick:
            tick()
        _time.sleep(TICK_S)

    fo.seek(0); fe.seek(0)
    out, err = fo.read(), fe.read()
    fo.close(); fe.close()

    res = {"bench": {}, "check": {}, "require": {}, "skip": {},
           "stderr": err.strip(), "timed_out": timed_out, "rc": rc,
           "error": None}
    for line in out.splitlines():
        m = LINE_RE.match(line.rstrip("\n"))
        if not m:
            continue
        tag, a, b = m.group(1), m.group(2).strip(), m.group(3).strip()
        if tag == "BENCH":
            try:
                res["bench"][a] = float(b)
            except ValueError:
                pass
        elif tag == "CHECK":
            res["check"][a] = b
        elif tag == "REQUIRE":
            res["require"][a] = b
        elif tag == "SKIP":
            res["skip"][a] = b
    return res


# --------------------------------------------------------------------------
# The value gate
# --------------------------------------------------------------------------

def as_float(s):
    try:
        return float(str(s).strip())
    except (TypeError, ValueError):
        return None


def norm_text(s):
    return re.sub(r"\s+", "", str(s)).rstrip(".").lower()


def checks_agree(values):
    """True when every recorded check value matches the first.

    Floats on relative error (the three systems do not round identically);
    everything else on normalised text.
    """
    vals = [v for v in values if v is not None and v != ""]
    if len(vals) < 2:
        return True                      # nothing to disagree with
    fs = [as_float(v) for v in vals]
    if all(f is not None for f in fs):
        ref = fs[0]
        scale = max(1.0, *(abs(f) for f in fs))
        return all(abs(f - ref) <= CHECK_RTOL * scale for f in fs)
    return len({norm_text(v) for v in vals}) == 1


# --------------------------------------------------------------------------
# Classification
# --------------------------------------------------------------------------

def classify(exp, runs):
    """One CASE record per label, in exactly one class."""
    rows = []          # cases
    labels = []
    for sysname in SYSTEMS:
        r = runs.get(sysname)
        if not r:
            continue
        # Check-only labels count too: a case that emitted a CHECK but no BENCH
        # is an INCOMPLETE finding, and dropping it would hide exactly the situation
        # where a head returns unevaluated.
        for lbl in list(r["bench"]) + list(r["skip"]) + list(r["check"]):
            if lbl not in labels:
                labels.append(lbl)

    mres = runs.get("mathilda") or {}
    for lbl in labels:
        times = {}
        for sysname in SYSTEMS:
            r = runs.get(sysname)
            if r and lbl in r["bench"]:
                times[sysname] = r["bench"][lbl]

        chk = {}
        for sysname in SYSTEMS:
            r = runs.get(sysname)
            if r and lbl in r["check"]:
                chk[sysname] = r["check"][lbl]

        row = {"experiment": exp["name"], "group": exp["group"], "label": lbl,
               "times": times, "checks": chk, "ratio": None, "note": ""}

        # 1. Absent -- declared by the .m side via benchIf/SKIP.
        if mres and lbl in mres.get("skip", {}):
            row["klass"] = "ABSENT"
            row["missing"] = mres["skip"][lbl]
            row["note"] = "Mathilda has no %s" % mres["skip"][lbl]
            rows.append(row)
            continue

        # 2. Timed out, or never produced a row at all.
        if "mathilda" not in times:
            row["klass"] = "INCOMPLETE"
            if mres.get("timed_out"):
                row["note"] = "timed out"
            elif mres.get("error"):
                row["note"] = mres["error"]
            else:
                row["note"] = "no row emitted (errored or returned unevaluated)"
            rows.append(row)
            continue

        # 3. The value gate outranks every timing.
        if not checks_agree(list(chk.values())):
            row["klass"] = "CHECK-FAIL"
            row["note"] = "check disagreement: " + "; ".join(
                "%s=%s" % (PRETTY[k], v) for k, v in sorted(chk.items()))
            row["times"] = {}            # DISCARD the timings, do not show them
            rows.append(row)
            continue

        # 4. Speed.
        others = {k: v for k, v in times.items() if k != "mathilda" and v > 0}
        if not others:
            row["klass"] = "AHEAD"
            row["note"] = "no baseline recorded"
            rows.append(row)
            continue
        best = min(others.values())
        mine = times["mathilda"]
        row["ratio"] = (mine / best) if best > 0 else None
        row["best_other"] = min(others, key=others.get)
        row["klass"] = "SLOWER" if (row["ratio"] and row["ratio"] >= SLOWER_AT) \
            else "AHEAD"
        rows.append(row)
    return rows


# --------------------------------------------------------------------------
# Host / version block -- the existing corpus's block is stale and misled a
# whole research pass, so this is recorded on every run.
# --------------------------------------------------------------------------

def host_block():
    def ver(mod):
        try:
            return __import__(mod).__version__
        except Exception:
            return "absent"

    mv = "unknown"
    try:
        mv = subprocess.run([MATHILDA, "-v"], capture_output=True, text=True,
                            timeout=60).stdout.strip() or "no output"
    except Exception as e:
        mv = "error: %s" % e

    pyver = "unknown"
    try:
        pyver = subprocess.run(
            [PYTHON, "-c",
             "import sys,numpy,scipy,sympy;"
             "print(sys.version.split()[0],numpy.__version__,"
             "scipy.__version__,sympy.__version__)"],
            capture_output=True, text=True, timeout=120).stdout.strip()
    except Exception as e:
        pyver = "error: %s" % e

    # THE BUILD'S OPTIONAL-DEPENDENCY STATE, and why it is recorded on every run.
    #
    # A missing optional dependency does not disable a feature here -- it
    # silently swaps in a slower algorithm.  With FFTW absent, `Fourier` falls
    # back to a naive O(n^2) DFT: a 65536-point transform measured at 4.88 s
    # against NumPy's ~1 ms.  Reported without this block, that row reads
    # "Mathilda is 5000x slower at FFT" and sends a week of work at
    # src/fourier.c, when the actual fix is `brew install fftw`.
    #
    # `Mathilda -v` cannot answer this: it lists what IS compiled in and says
    # nothing about what is missing, which is precisely how this stayed
    # invisible. The makefile's own detection warnings are the exact answer.
    missing = []
    try:
        mk = subprocess.run(["make", "-n"], cwd=ROOT, capture_output=True,
                            text=True, timeout=180)
        for line in (mk.stdout + mk.stderr).splitlines():
            if "not detected" in line:
                missing.append(line.split(":", 2)[-1].strip())
    except Exception as e:
        missing.append("could not probe: %s" % e)

    return {
        "date": datetime.date.today().isoformat(),
        "missing_optional_deps": missing,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor() or platform.machine(),
        "mathilda": mv,
        "mathilda_bin": MATHILDA,
        "python_bin": PYTHON,
        "python_libs": pyver,
        "wolframscript": WOLFRAM or "absent",
        "runner_python": sys.version.split()[0],
        "local_versions": {m: ver(m) for m in
                           ("numpy", "scipy", "sympy", "mpmath", "numba",
                            "networkx")},
    }


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------

# Below this, a baseline timing is at the clock's resolution and a ratio
# computed from it carries no real precision. 0.005 ms is where the corpus's
# reported values start rounding to 0.00x.
RATIO_FLOOR_MS = 0.005


def fmt_ratio(r, best_ms):
    """Render a ratio, refusing to imply precision the baseline cannot support.

    A case where the baseline ran in 0.001 ms is not "6097090.0x behind" to seven
    significant figures -- 0.001 ms is the floor of this formatter and near
    AbsoluteTiming's granularity, so the true figure could be several times
    either side. Those render as ">N x" with a footnote. The finding is unharmed
    (6.10 s for ConnectedComponents is indefensible at any precision); what is
    removed is the false precision.
    """
    if r is None:
        return "—"
    if best_ms is not None and best_ms < RATIO_FLOOR_MS:
        # Round down to one significant figure and mark it as a lower bound.
        import math
        if r <= 0:
            return "—"
        mag = 10 ** int(math.floor(math.log10(r)))
        return "**>%s×**\u2009*" % ("{:,}".format(int(mag)))
    if r >= 100:
        return "**%.0f×**" % r
    return "**%.1f×**" % r


def median(xs):
    xs = sorted(x for x in xs if x is not None)
    if not xs:
        return None
    m = len(xs) // 2
    return xs[m] if len(xs) % 2 else 0.5 * (xs[m - 1] + xs[m])


def ms(v):
    if v is None:
        return "—"
    if v >= 1000:
        return "%.2f s" % (v / 1000.0)
    if v >= 1:
        return "%.1f ms" % v
    return "%.3f ms" % v


def coverage(all_requires):
    """present / total over every head the .m files declared via require[]."""
    present = sorted(h for h, s in all_requires.items() if s == "present")
    absent = sorted(h for h, s in all_requires.items() if s != "present")
    total = len(present) + len(absent)
    pct = (100.0 * len(present) / total) if total else 0.0
    return present, absent, pct



def load_prev_run(today):
    """The most recent FULL run before today, for week-over-week comparison.

    Partial runs are skipped: a `--only` subset would make the trend read as if
    two dozen cases had vanished.
    """
    outdir = os.path.join(BENCH, "results")
    if not os.path.isdir(outdir):
        return None
    cur = os.path.join(outdir, "%s.json" % today)
    cands = [os.path.join(outdir, n) for n in os.listdir(outdir)
             if n.endswith(".json") and "partial" not in n]
    # By MODIFICATION TIME, not filename: "2026-08-04.1.json" sorts BEFORE
    # "2026-08-04.json" lexically ('.' < 'j'), so name order picks the wrong one.
    # Excluding today's canonical file means a same-day re-run compares against
    # the archived earlier run rather than against itself.
    for path in sorted(cands, key=lambda q: -os.path.getmtime(q)):
        if os.path.abspath(path) == os.path.abspath(cur):
            continue
        try:
            with open(path) as f:
                d = json.load(f)
            if d.get("rows"):
                return d
        except Exception:
            continue
    return None


def trend_section(rows, requires, prev):
    """Week-over-week: what got better, what got worse, what is new."""
    if not prev:
        return ["## Trend", "",
                "_No previous full run to compare against. The next run will diff "
                "against this one — that is what makes this a weekly job rather "
                "than a snapshot._", ""]

    def key(r):
        return (r["experiment"], r["label"])

    now = {key(r): r for r in rows}
    was = {key(r): r for r in prev["rows"]}
    pdate = prev["host"].get("date", "?")

    improved, regressed, newly, gone = [], [], [], []
    for k, r in now.items():
        o = was.get(k)
        if o is None:
            newly.append(r)
            continue
        a, b = o.get("ratio"), r.get("ratio")
        if a and b and a > 0:
            change = b / a
            if change <= 0.8:
                improved.append((r, a, b))
            elif change >= 1.25:
                regressed.append((r, a, b))
        elif o["klass"] != r["klass"]:
            if r["klass"] in ("AHEAD",) and o["klass"] != "AHEAD":
                improved.append((r, None, None))
            elif o["klass"] == "AHEAD" and r["klass"] != "AHEAD":
                regressed.append((r, None, None))
    for k, o in was.items():
        if k not in now:
            gone.append(o)

    pc = {kl: sum(1 for r in prev["rows"] if r["klass"] == kl)
          for kl in ("SLOWER", "AHEAD", "ABSENT", "INCOMPLETE", "CHECK-FAIL")}
    nc = {kl: sum(1 for r in rows if r["klass"] == kl)
          for kl in ("SLOWER", "AHEAD", "ABSENT", "INCOMPLETE", "CHECK-FAIL")}
    _, pabs, ppct = coverage(prev.get("requires", {}))
    _, nabs, npct = coverage(requires)

    L = ["## Trend since %s" % pdate, ""]
    L.append("| | then | now | Δ |")
    L.append("|---|---|---|---|")
    L.append("| coverage | %.1f%% | %.1f%% | %+.1f pp |" % (ppct, npct, npct - ppct))
    for kl in ("SLOWER", "AHEAD", "ABSENT", "INCOMPLETE", "CHECK-FAIL"):
        L.append("| %s | %d | %d | %+d |" % (kl, pc[kl], nc[kl], nc[kl] - pc[kl]))
    L.append("")

    def _tbl(title, items, note):
        out = ["### %s (%d)" % (title, len(items)), ""]
        if not items:
            out += ["_None._", ""]
            return out
        out.append(note)
        out.append("")
        out.append("| case | then | now |")
        out.append("|---|---|---|")
        for r, a, b in sorted(items,
                              key=lambda t: -(t[2] / t[1]) if t[1] and t[2] else 0)[:12]:
            out.append("| `%s` · %s | %s | %s |" % (
                r["experiment"], r["label"],
                ("%.1f×" % a) if a else "—", ("%.1f×" % b) if b else "—"))
        out.append("")
        return out

    L += _tbl("Regressed", regressed,
              "At least 25% worse relative to the best baseline. **This is the "
              "half of the diff worth reading first** — a benchmark job that only "
              "reports wins is a marketing document.")
    L += _tbl("Improved", improved, "At least 20% better relative to the best baseline.")
    if newly or gone:
        L.append("### Coverage of the suite itself")
        L.append("")
        if newly:
            L.append("- %d new case(s) since %s." % (len(newly), pdate))
        if gone:
            L.append("- %d case(s) present then and missing now — check whether an "
                     "experiment stopped emitting rather than assuming it passed."
                     % len(gone))
        L.append("")
    return L


def history_record(rows, requires, host, elapsed):
    """The machine-readable summary of one full run.

    Kept as data, not prose. `history.jsonl` is the source of truth — one JSON
    object per line, append-only, so it diffs cleanly in git, can be read by
    `json.loads` line by line without parsing markdown, and can be plotted
    directly. `HISTORY.md` is a *rendered view* of it, regenerated every run.

    Per-group medians are stored per baseline so an area's trend can be plotted
    on its own: "is group C closing on NumPy" is a question about one series
    here, not something to be reconstructed by re-parsing old reports.
    """
    present, absent, pct = coverage(requires)
    counts = {kl: sum(1 for r in rows if r["klass"] == kl)
              for kl in ("SLOWER", "AHEAD", "ABSENT", "INCOMPLETE", "CHECK-FAIL")}

    def rel(subset, sysname):
        v = [r["times"]["mathilda"] / r["times"][sysname] for r in subset
             if "mathilda" in r["times"] and r["times"].get(sysname, 0) > 0]
        return median(v)

    groups = {}
    for g in sorted({r["group"] for r in rows}):
        gr = [r for r in rows if r["group"] == g]
        groups[g] = {
            "cases": len(gr),
            "slower": sum(1 for r in gr if r["klass"] == "SLOWER"),
            "ahead": sum(1 for r in gr if r["klass"] == "AHEAD"),
            "absent": sum(1 for r in gr if r["klass"] == "ABSENT"),
            "incomplete": sum(1 for r in gr if r["klass"] in
                              ("INCOMPLETE", "CHECK-FAIL")),
            "median_vs_wolfram": rel(gr, "wolfram"),
            "median_vs_python": rel(gr, "python"),
        }

    commit = "unknown"
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                                capture_output=True, text=True,
                                timeout=30).stdout.strip() or "unknown"
    except Exception:
        pass

    return {
        "date": host["date"],
        "commit": commit,
        "cases": len(rows),
        "coverage_pct": round(pct, 2),
        "heads_present": len(present),
        "heads_absent": len(absent),
        "counts": counts,
        "median_vs_wolfram": rel(rows, "wolfram"),
        "median_vs_python": rel(rows, "python"),
        "groups": groups,
        "wall_min": round(elapsed / 60.0, 2),
        "missing_optional_deps": host.get("missing_optional_deps", []),
        "mathilda": host.get("mathilda", ""),
        "platform": host.get("platform", ""),
    }


def append_history(rows, requires, host, elapsed):
    """Append this run to history.jsonl, then re-render HISTORY.md from it."""
    jsonl = os.path.join(BENCH, "history.jsonl")
    rec = history_record(rows, requires, host, elapsed)

    # A same-day re-run replaces its own row rather than adding a second.
    kept = []
    if os.path.exists(jsonl):
        with open(jsonl) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    if json.loads(line).get("date") != rec["date"]:
                        kept.append(line)
                except Exception:
                    kept.append(line)
    kept.append(json.dumps(rec, sort_keys=True))
    with open(jsonl, "w") as f:
        f.write("\n".join(kept) + "\n")

    render_history_md([json.loads(l) for l in kept])


def render_history_md(recs):
    """HISTORY.md is a VIEW of history.jsonl — never edited, always regenerated."""
    L = ["# Benchmark history", "",
         "Rendered from [`history.jsonl`](history.jsonl), which is the source of "
         "truth — one JSON object per full run, append-only. **Do not edit this "
         "file**; it is regenerated by `make bench-gap`. Query or plot the JSONL "
         "instead:", "",
         "```bash",
         "python3 -c \"import json;[print(json.loads(l)['date'], "
         "json.loads(l)['median_vs_python']) for l in open('benchmarks/history.jsonl')]\"",
         "```", "",
         "Medians, not means — one 600000× case would make a mean meaningless.",
         "",
         "| date | commit | cases | coverage | slower | ahead | absent | "
         "incomplete | vs Mathematica | vs Python | wall |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]

    def r2(v):
        return ("%.2f×" % v) if isinstance(v, (int, float)) else "—"

    for rec in recs:
        c = rec.get("counts", {})
        L.append("| %s | `%s` | %d | %.1f%% | %d | %d | %d | %d | %s | %s | %.1f min |"
                 % (rec.get("date", "?"), rec.get("commit", "?"),
                    rec.get("cases", 0), rec.get("coverage_pct", 0.0),
                    c.get("SLOWER", 0), c.get("AHEAD", 0), c.get("ABSENT", 0),
                    c.get("INCOMPLETE", 0) + c.get("CHECK-FAIL", 0),
                    r2(rec.get("median_vs_wolfram")),
                    r2(rec.get("median_vs_python")), rec.get("wall_min", 0.0)))
    L.append("")

    if recs:
        L.append("## By area, over time")
        L.append("")
        L.append("Median ratio per area. The question this answers is which "
                 "subsystem is actually closing, as opposed to which one got "
                 "attention.")
        L.append("")
        areas = sorted({g for rec in recs for g in rec.get("groups", {})})
        L.append("| date | " + " | ".join("%s<br>(vs W / vs Py)" % a for a in areas) + " |")
        L.append("|---" * (len(areas) + 1) + "|")
        for rec in recs:
            cells = []
            for a in areas:
                g = rec.get("groups", {}).get(a, {})
                cells.append("%s / %s" % (r2(g.get("median_vs_wolfram")),
                                          r2(g.get("median_vs_python"))))
            L.append("| %s | %s |" % (rec.get("date", "?"), " | ".join(cells)))
        L.append("")

    with open(os.path.join(BENCH, "HISTORY.md"), "w") as f:
        f.write("\n".join(L) + "\n")


def write_report(path, rows, requires, host, exps, elapsed):
    present, absent, pct = coverage(requires)
    by = {k: [r for r in rows if r["klass"] == k]
          for k in ("SLOWER", "AHEAD", "ABSENT", "INCOMPLETE", "CHECK-FAIL")}
    slower = sorted(by["SLOWER"], key=lambda r: -(r["ratio"] or 0))

    L = []
    L.append("# Benchmark gap report — %s" % host["date"])
    L.append("")
    L.append("Generated by `make bench-gap` (`benchmarks/run_all.py`). Do not edit "
             "by hand.")
    L.append("")
    L.append("## Headline")
    L.append("")
    L.append("| | |")
    L.append("|---|---|")
    L.append("| Experiments run | %d |" % len(exps))
    L.append("| Cases measured | %d |" % len(rows))
    L.append("| **Function coverage** | **%.1f%%** (%d of %d declared heads present) |"
             % (pct, len(present), len(present) + len(absent)))
    L.append("| Slower than best baseline | %d |" % len(by["SLOWER"]))
    L.append("| Ahead or within %.1f× | %d |" % (SLOWER_AT, len(by["AHEAD"])))
    L.append("| Absent (feature work) | %d |" % len(by["ABSENT"]))
    L.append("| Incomplete / timed out | %d |" % len(by["INCOMPLETE"]))
    L.append("| Check disagreements | %d |" % len(by["CHECK-FAIL"]))
    L.append("| Wall clock | %.1f min |" % (elapsed / 60.0))
    L.append("")
    L.append("`SLOWER` cases carry a ratio and are kernel work. `ABSENT` cases carry "
             "**no ratio** and are feature work — see `ABSENT.md`. The two are never "
             "pooled: averaging \"we do not have `DSolve`\" into a speed number "
             "would turn a missing algorithm into a fake 40×.")
    L.append("")

    L += trend_section(rows, requires, load_prev_run(host["date"]))

    L.append("## Ranked gaps — where the week's work is")
    L.append("")
    if not slower:
        L.append("_No case is %.1f× or more behind its best baseline._" % SLOWER_AT)
    else:
        L.append("| × behind | case | Mathilda | Mathematica | Python | vs |")
        L.append("|---|---|---|---|---|---|")
        any_floored = False
        for r in slower:
            t = r["times"]
            others = [v for k, v in t.items() if k != "mathilda"]
            bo = min(others) if others else None
            if bo is not None and bo < RATIO_FLOOR_MS:
                any_floored = True
            L.append("| %s | `%s` · %s | %s | %s | %s | %s |" % (
                fmt_ratio(r["ratio"], bo), r["experiment"], r["label"],
                ms(t.get("mathilda")), ms(t.get("wolfram")), ms(t.get("python")),
                PRETTY.get(r.get("best_other"), "—")))
        if any_floored:
            L.append("")
            L.append("\u2009* baseline ran under %.3f ms — at clock resolution, so "
                     "the ratio is a lower bound, not a measurement. The gap is "
                     "real; its magnitude is not precise." % RATIO_FLOOR_MS)
    L.append("")

    L.append("## Where Mathilda leads or ties")
    L.append("")
    ahead = sorted(by["AHEAD"], key=lambda r: (r["ratio"] or 0))
    if not ahead:
        L.append("_None._")
    else:
        L.append("| case | Mathilda | Mathematica | Python | ratio |")
        L.append("|---|---|---|---|---|")
        for r in ahead:
            t = r["times"]
            L.append("| `%s` · %s | %s | %s | %s | %s |" % (
                r["experiment"], r["label"], ms(t.get("mathilda")),
                ms(t.get("wolfram")), ms(t.get("python")),
                ("%.2f×" % r["ratio"]) if r["ratio"] else "—"))
    L.append("")

    if by["CHECK-FAIL"]:
        L.append("## Check disagreements — timings withheld")
        L.append("")
        L.append("These cases are **not** timing comparisons: the systems did not "
                 "agree on the answer, so any number here would be comparing "
                 "different programs. Fix the disagreement, then re-run.")
        L.append("")
        L.append("| case | detail |")
        L.append("|---|---|")
        for r in by["CHECK-FAIL"]:
            L.append("| `%s` · %s | %s |" % (r["experiment"], r["label"],
                                             r["note"]))
        L.append("")

    L.append("## Where the slowness is, by area")
    L.append("")
    L.append("Median ratios, not means: one 600000× case would drag a mean so far "
             "that the column would say nothing about the other sixty. A median "
             "of 1.0 means the typical case in that area is at parity.")
    L.append("")
    L.append("| area | cases | slower | ahead | absent | incomplete | "
             "median vs Mathematica | median vs Python |")
    L.append("|---|---|---|---|---|---|---|---|")

    def _rel(gr, sysname):
        """Median of (Mathilda / that system) over cases where both timed."""
        out = []
        for r in gr:
            t = r["times"]
            if "mathilda" in t and t.get(sysname, 0) > 0:
                out.append(t["mathilda"] / t[sysname])
        return median(out)

    for g in sorted({r["group"] for r in rows}):
        gr = [r for r in rows if r["group"] == g]
        mw, mp = _rel(gr, "wolfram"), _rel(gr, "python")
        L.append("| %s | %d | %d | %d | %d | %d | %s | %s |" % (
            g, len(gr),
            sum(1 for r in gr if r["klass"] == "SLOWER"),
            sum(1 for r in gr if r["klass"] == "AHEAD"),
            sum(1 for r in gr if r["klass"] == "ABSENT"),
            sum(1 for r in gr if r["klass"] in ("INCOMPLETE", "CHECK-FAIL")),
            ("%.2f×" % mw) if mw else "—",
            ("%.2f×" % mp) if mp else "—"))
    L.append("")
    L.append("Read the two right-hand columns against each other. Behind "
             "Mathematica is a gap against a *competitor* — the same evaluation "
             "semantics, the same automatic compilation — so it is an engineering "
             "gap. Behind Python is a gap against *the machine*: on the dense "
             "linear-algebra areas all three call the same Accelerate BLAS, so "
             "any spread there is pure overhead. A case can read acceptably "
             "against one and badly against the other, and that difference is "
             "usually the finding.")
    L.append("")

    L.append("### The worst area, by count")
    L.append("")
    worst = sorted(({"g": g,
                     "n": sum(1 for r in rows
                              if r["group"] == g and r["klass"] == "SLOWER")}
                    for g in {r["group"] for r in rows}),
                   key=lambda d: -d["n"])
    for d in worst[:2]:
        gr = [r for r in rows if r["group"] == d["g"] and r["klass"] == "SLOWER"]
        top = sorted(gr, key=lambda r: -(r["ratio"] or 0))[:3]
        L.append("- **%s** — %d slower cases. Worst: %s" % (
            d["g"], d["n"],
            "; ".join("`%s` (%s)" % (r["label"][:38],
                                     ("%.0f×" % r["ratio"]) if r["ratio"] else "?")
                      for r in top)))
    L.append("")

    L.append("## Method")
    L.append("")
    L.append("- One untimed warm-up, then the **minimum** of 3 timed runs "
             "(`AbsoluteTiming` / `perf_counter`; never `Timing[]`, which sums CPU "
             "time over threads). Cases already expensive on their own use one "
             "timed run and no warm-up.")
    L.append("- **One process per experiment per system**, cwd set to the "
             "experiment folder.")
    L.append("- The same `.m` source runs in Mathilda and Mathematica. The `.py` is "
             "a separate implementation of the same algorithm in the same order.")
    L.append("- **Values are checked, not just timed.** A case whose systems "
             "disagree is reported `CHECK-FAIL` with its timing discarded.")
    L.append("- Exit codes are ignored: Mathilda exits 0 after errors, so "
             "classification comes from parsed cases plus stderr.")
    L.append("")
    L.append("### Caveat that applies to group A only")
    L.append("")
    L.append("Symbolic timings are frequently dominated by **algorithm choice, not "
             "execution speed**. A group-A case reading 40× behind sympy usually "
             "means *sympy has a method we do not*, which is different work from "
             "the numeric findings (nearly all of which were \"the fast path "
             "existed and was not taken\"). Read group A cases together with "
             "`ABSENT.md`.")
    L.append("")

    if host.get("missing_optional_deps"):
        L.append("## ⚠ Build warning — read before trusting any row")
        L.append("")
        L.append("This build is missing optional dependencies. A missing optional "
                 "dependency here does not disable a feature — it **silently "
                 "substitutes a slower algorithm**, so the affected cases measure "
                 "the fallback, not Mathilda's real implementation. Do not open "
                 "kernel work on these cases until the dependency is installed and "
                 "the tree rebuilt.")
        L.append("")
        for w in host["missing_optional_deps"]:
            L.append("- %s" % w)
        L.append("")

    L.append("## Host")
    L.append("")
    L.append("| | |")
    L.append("|---|---|")
    for k in ("platform", "machine", "processor", "mathilda", "python_bin",
              "python_libs", "wolframscript", "runner_python"):
        L.append("| %s | `%s` |" % (k, host[k]))
    L.append("")

    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


def write_absent(path, rows, requires, host):
    present, absent, pct = coverage(requires)
    absent_rows = [r for r in rows if r["klass"] == "ABSENT"]
    incomplete = [r for r in rows if r["klass"] == "INCOMPLETE"]

    L = []
    L.append("# Absence register — %s" % host["date"])
    L.append("")
    L.append("Generated by `make bench-gap`. Do not edit by hand.")
    L.append("")
    L.append("What Mathilda **does not have**, kept strictly apart from what "
             "Mathilda has and runs slowly. Nothing here carries a speed ratio: a "
             "missing function is feature work, and giving it a ratio would launder "
             "a capability gap into a performance number.")
    L.append("")
    L.append("**Coverage: %.1f%% — %d of %d declared heads present, %d absent.**"
             % (pct, len(present), len(present) + len(absent), len(absent)))
    L.append("")

    L.append("## Absent heads")
    L.append("")
    if not absent:
        L.append("_None._")
    else:
        L.append("| head | blocks |")
        L.append("|---|---|")
        for h in absent:
            blockers = sorted({r["label"] for r in absent_rows
                               if r.get("missing") == h})
            L.append("| `%s` | %s |" % (
                h, ", ".join(blockers) if blockers else "_declared, not yet a row_"))
    L.append("")

    if absent_rows:
        L.append("## Cases that could not run")
        L.append("")
        L.append("| experiment | case | missing |")
        L.append("|---|---|---|")
        for r in sorted(absent_rows, key=lambda r: (r["experiment"], r["label"])):
            L.append("| `%s` | %s | `%s` |" % (r["experiment"], r["label"],
                                               r.get("missing", "?")))
        L.append("")

    if incomplete:
        L.append("## Incomplete — the head exists but the case did not finish")
        L.append("")
        L.append("Distinct from absence: these are correctness or coverage bugs, "
                 "not missing features. A head that returns **unevaluated** is the "
                 "dangerous case, because it costs almost no time and so reads as a "
                 "win unless something checks the value.")
        L.append("")
        L.append("| experiment | case | why |")
        L.append("|---|---|---|")
        for r in sorted(incomplete, key=lambda r: (r["experiment"], r["label"])):
            L.append("| `%s` | %s | %s |" % (r["experiment"], r["label"],
                                             r["note"]))
        L.append("")

    L.append("## Present heads (%d)" % len(present))
    L.append("")
    L.append(", ".join("`%s`" % h for h in present) if present else "_None._")
    L.append("")

    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    import time as _t

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", help="comma-separated experiment numbers or slugs")
    ap.add_argument("--system", default="mathilda,wolfram,python",
                    help="comma-separated subset of mathilda,wolfram,python")
    ap.add_argument("--timeout", type=int, default=300,
                    help="per-file wall-clock seconds (default 300)")
    ap.add_argument("--list", action="store_true", help="list experiments and exit")
    ap.add_argument("--check-labels", action="store_true",
                    help="fail if an experiment's .m and .py emit different labels")
    ap.add_argument("--json", help="raw results path (default results/<date>.json)")
    ap.add_argument("-v", "--verbose", action="store_true",
                    help="print a line per experiment and per system instead of "
                         "a single-line progress bar")
    ap.add_argument("--from-json", metavar="PATH",
                    help="re-render REPORT.md/ABSENT.md from a saved results "
                         "file and exit, without re-measuring anything. Use "
                         "after changing report wording or thresholds.")
    args = ap.parse_args()

    if args.from_json:
        with open(args.from_json) as f:
            d = json.load(f)
        rows, reqs, host = d["rows"], d["requires"], d["host"]
        exps_seen = sorted({r["experiment"] for r in rows})
        write_report(os.path.join(BENCH, "REPORT.md"), rows, reqs, host,
                     exps_seen, d.get("elapsed_s", 0.0))
        write_absent(os.path.join(BENCH, "ABSENT.md"), rows, reqs, host)
        print("re-rendered from %s (no history row appended): "
              "%d cases, %d experiments"
              % (args.from_json, len(rows), len(exps_seen)), file=sys.stderr)
        return 0

    only = set(x.strip() for x in args.only.split(",")) if args.only else None
    want = [s.strip() for s in args.system.split(",") if s.strip()]
    exps = discover(only)

    if args.list:
        for e in exps:
            print("%s  %-34s %s" % (e["num"], e["slug"], e["group"]))
        return 0
    if not exps:
        print("no experiments found under %s" % BENCH, file=sys.stderr)
        return 2

    if "wolfram" in want and not WOLFRAM:
        print("note: wolframscript not found — the Mathematica column will read "
              "'—'", file=sys.stderr)

    verbose = args.verbose
    prev = load_prev_durations()
    total = len(exps)
    print("%d experiments x %s  (timeout %ds/file)"
          % (total, "+".join(want), args.timeout), file=sys.stderr)
    if prev:
        est = sum(prev.get(e["name"], 0) for e in exps)
        print("last run's total for these: %s — that is the ETA basis"
              % fmt_dur(est), file=sys.stderr)
    else:
        print("no previous run to calibrate against; the ETA will be rough "
              "until a few experiments have finished", file=sys.stderr)
    print("", file=sys.stderr)

    t_start = _t.time()
    host = host_block()
    all_rows, all_requires, label_problems, raw = [], {}, [], {}
    durations, done_names = {}, []   # done_names also feeds eta()

    bar_ui = Progress(total, exps, prev, verbose)
    bar_ui.cap = args.timeout * len(want)   # a hung experiment's ceiling
    bar_ui.start(t_start)

    for i, e in enumerate(exps, 1):
        t_exp = _t.time()
        bar_ui.begin(i, e["slug"])
        if verbose:
            bar_ui.note("[%2d/%d] %s" % (i, total, e["slug"]))
        runs = {}
        for sysname in SYSTEMS:
            if sysname not in want:
                continue
            t_sys = _t.time()
            bar_ui.system(sysname)
            r = run_one(sysname, e, args.timeout, tick=bar_ui.draw)
            if r is None:
                continue
            runs[sysname] = r
            n = len(r["bench"])
            extra = " TIMEOUT" if r["timed_out"] else ""
            if verbose or r["timed_out"]:
                bar_ui.note("        %-11s %3d cases %7s%s"
                            % (sysname, n, fmt_dur(_t.time() - t_sys), extra))
        durations[e["name"]] = _t.time() - t_exp
        done_names.append(e["name"])
        bar_ui.finish_experiment(e["name"])

        # Mathilda's require[] lines are what the coverage percentage counts:
        # the question is what MATHILDA has, not what sympy has.
        mr = runs.get("mathilda")
        if mr:
            for h, state in mr["require"].items():
                if all_requires.get(h) != "absent":
                    all_requires[h] = state

        # Label-set agreement between the two halves of the experiment.
        if "mathilda" in runs and "python" in runs:
            a = set(runs["mathilda"]["bench"]) | set(runs["mathilda"]["skip"])
            b = set(runs["python"]["bench"]) | set(runs["python"]["skip"])
            if a != b:
                label_problems.append((e["name"], sorted(a - b), sorted(b - a)))

        all_rows.extend(classify(e, runs))
        raw[e["name"]] = {s: runs[s] for s in runs}

    elapsed = _t.time() - t_start
    bar_ui.close()

    # A PARTIAL RUN MUST NOT CLOBBER THE WEEKLY REPORT.  `--only` / `--system`
    # produce a subset, and the first time this was used for a 3-experiment
    # smoke test it overwrote the full suite's REPORT.md and its results JSON --
    # destroying both the committed weekly artefact and the per-experiment
    # timings the ETA calibrates against.  Subsets therefore write `.partial`
    # names and leave the canonical files alone.
    partial = bool(only) or set(want) != set(SYSTEMS)
    suffix = "-partial" if partial else ""
    outdir = os.path.join(BENCH, "results")
    os.makedirs(outdir, exist_ok=True)
    jpath = args.json or os.path.join(outdir, "%s%s.json" % (host["date"], suffix))

    # A second full run on the same day used to OVERWRITE the first, destroying
    # the very run the trend should diff against — so re-running to check a fix
    # silently cost you the before-picture. Rotate the old one aside instead.
    # `<date>.json` stays the canonical latest; `<date>.1.json` and up are the
    # earlier runs of that day, newest-highest.
    if not partial and not args.json and os.path.exists(jpath):
        k = 1
        while os.path.exists(os.path.join(outdir, "%s.%d.json" % (host["date"], k))):
            k += 1
        os.rename(jpath, os.path.join(outdir, "%s.%d.json" % (host["date"], k)))
        print("  previous run of today archived as results/%s.%d.json"
              % (host["date"], k), file=sys.stderr)
    with open(jpath, "w") as f:
        json.dump({"host": host, "rows": all_rows, "requires": all_requires,
                   "raw": raw, "elapsed_s": elapsed,
                   "durations": durations, "partial": partial,
                   "label_problems": label_problems}, f, indent=1)

    rname = "REPORT.partial.md" if partial else "REPORT.md"
    aname = "ABSENT.partial.md" if partial else "ABSENT.md"
    write_report(os.path.join(BENCH, rname), all_rows, all_requires, host,
                 exps, elapsed)
    write_absent(os.path.join(BENCH, aname), all_rows, all_requires, host)
    if not partial:
        # Only a full run earns a history row; a subset would distort the arc.
        append_history(all_rows, all_requires, host, elapsed)

    present, absent, pct = coverage(all_requires)
    counts = {k: sum(1 for r in all_rows if r["klass"] == k)
              for k in ("SLOWER", "AHEAD", "ABSENT", "INCOMPLETE", "CHECK-FAIL")}
    print("\n%d cases in %.1f min — coverage %.1f%% (%d/%d heads)"
          % (len(all_rows), elapsed / 60.0, pct, len(present),
             len(present) + len(absent)), file=sys.stderr)
    print("  " + "  ".join("%s=%d" % (k, v) for k, v in counts.items()),
          file=sys.stderr)
    slowest = sorted(durations.items(), key=lambda kv: -kv[1])[:5]
    if slowest:
        print("  slowest: " + ", ".join("%s %s" % (n, fmt_dur(d))
                                        for n, d in slowest),
              file=sys.stderr)
    print("  wrote benchmarks/%s, benchmarks/%s, %s%s"
          % (rname, aname, os.path.relpath(jpath, ROOT),
             "   (PARTIAL run: canonical REPORT.md left untouched)"
             if partial else ""), file=sys.stderr)

    if label_problems:
        print("\nLABEL MISMATCH between the .m and .py halves:", file=sys.stderr)
        for name, onlym, onlypy in label_problems:
            print("  %s" % name, file=sys.stderr)
            for l in onlym:
                print("      only .m : %s" % l, file=sys.stderr)
            for l in onlypy:
                print("      only .py: %s" % l, file=sys.stderr)
        if args.check_labels:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
