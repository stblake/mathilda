#!/usr/bin/env python3
"""
nd_fastpath_sweep.py -- does every builtin that CAN take a machine array
actually reach the buffer, measured rather than declared?

WHY A FOURTH TOOL. Three audits already guard the packed surface and all three
were green on the day `Commonest` was found costing 880 ms where `Tally` of the
identical 10^7 buffer cost 21.5 ms:

    check_packed_aware.py    reads the source for an NDArray dispatch and diffs
                             it against pack.c's AWARE list.  Commonest had NO
                             dispatch, so there was nothing to diff -- a head
                             with no fast path at all is invisible to it by
                             construction.
    check_array_exactness.py looks at the element HEADS of a result.  Commonest
                             answered with exactly the right heads, slowly.
    nd_surface_audit.py      requires List / packed / NDArray to agree.  All
                             three agreed, equally slowly.
    numeric_sweep.py         measures, and is the right shape -- but its probe
                             list is hand-written, and nobody had written a
                             `commonest` probe.

Every one of those either reads the source, or reads a curated list of
expressions.  The class they share a blind spot for is "a head nobody thought
to name".  So this tool names nobody: it enumerates every builtin the system
registers, discovers by TRIAL which call shapes that head accepts, and then
measures the shapes that work.  A head is a finding because of what the clock
says, not because it appears in a table someone maintained.

TWO PASSES, and the fast one is the detector.

  --gate-only   THE DETECTOR, in three phases.  MATHILDA_PACK_DIAG=gate makes
                pack.c count, per head, how many elements the transparency gate
                MATERIALISED -- one boxed Expr per element, built before the
                builtin ever ran.  That is exactly the mechanism a missing fast
                path shows up as, and it is a COUNT, not a duration:
                deterministic, immune to whatever else is on the machine, and
                complete in minutes where timing every head twice takes hours.

                  1. DISCOVER, at n = 8, which shapes each head accepts.
                  2. GATE only those shapes.  A head handed a shape it does not
                     take trips the post-gate too, and that is not a finding.
                  3. VERIFY: re-ask every survivor in ONE process, one explicit
                     expression per head.  The counter is per PROCESS and phase
                     2 batches 16 probes into one, so a head can be credited
                     with a neighbour's materialisation -- Exp, Gamma, Variance
                     and Clip were all named by phase 2 and none of them
                     materialises when asked on its own.

                Exits non-zero on a head that materialises and is not in the
                checked-in OFF_BUFFER list, so the standing backlog is reported
                and a NEW regression fails.

  (default)     THE SEVERITY.  Times the discovered shapes, so the findings can
                be ranked and a fix can be measured.  Run it over the heads the
                gate pass named, on an idle machine -- a timing tool sharing a
                machine with anything else is measuring the other thing too.

THE MEASUREMENT.  Two numbers per (head, shape), both wall clock via
AbsoluteTiming (never Timing[], which sums CPU time across threads):

    ns/elem   time divided by the element count of the argument.  A head on the
              buffer costs single-digit ns/elem; a head boxing one Expr per
              element costs hundreds.
    nopack    the same expression under MATHILDA_NO_PACK=1, which turns the
              packed representation off process-wide, over the ratio
              t_nopack / t_packed.

Neither alone is enough, which is the whole reason both are here:

  * ns/elem alone cannot tell a missing fast path from an expensive function.
    Zeta over 10^6 doubles is slow because Zeta is slow.
  * the ratio alone cannot tell a missing fast path from a head that never
    touches elements.  Length is ratio 1.0 and costs nothing either way.

Together they separate cleanly.  A head that is BOTH expensive per element AND
indifferent to whether its input was packed is not using a buffer -- there is
nothing else that shape of result can mean.  That is exactly the signature
Commonest had, and the signature Length, Zeta and Sin all fail to have.

Usage:
    python3 tools/nd_fastpath_sweep.py --gate-only       # the detector
    python3 tools/nd_fastpath_sweep.py                   # discover, measure, rank
    python3 tools/nd_fastpath_sweep.py --only Commonest,Tally
    python3 tools/nd_fastpath_sweep.py --discover-only   # which shapes evaluate
    python3 tools/nd_fastpath_sweep.py --n 1000000
    python3 tools/nd_fastpath_sweep.py --json out.json
    python3 tools/nd_fastpath_sweep.py --all             # print every row
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
MATHILDA = os.environ.get("MATHILDA_BIN", os.path.join(ROOT, "Mathilda"))

REPS = 2
# Probes per process.  Discovery runs at n = 8 and is nearly hang-free, so it
# batches hard: the whole cost of that phase is process spawns (7656 probes at
# 24 per batch was 319 of them, and ~25 of the 35 minutes the first full run
# took).  A batch that does time out is BISECTED rather than retried one probe
# at a time, so a big batch is not a big risk -- finding one bad probe in 128
# costs 8 re-runs, not 128.
DISCOVER_BATCH = 128
TIME_BATCH = 6
# Per-batch ceilings.  Both are deliberately tight: a probe that cannot answer
# in this long is itself the finding ("does not complete"), and waiting longer
# does not make it a better one.  The first full run used 90 s and spent most of
# its wall clock in single-probe retries behind pathological shapes such as
# CountsBy[vi, wi] -- two 200000-element integer vectors where the second is a
# key list, which is quadratic and was never going to answer.
TIMEOUT = int(os.environ.get("SWEEP_TIMEOUT", "25"))
DISCOVER_TIMEOUT = 10

# Size used for the SKIP_EXPLOSIVE heads.  ns/elem is scale-free, so a smaller
# vector still answers "is this on the buffer"; it only widens the error bar.
SMALL_N = 5000
SMALL_MDIM = 70


# ---------------------------------------------------------------------------
# The head inventory: every name the system registers, from the source.
# ---------------------------------------------------------------------------

def _walk_c():
    for dirpath, dirnames, filenames in os.walk(SRC):
        dirnames[:] = [d for d in dirnames if d != "external"]
        for fn in filenames:
            if fn.endswith(".c"):
                yield os.path.join(dirpath, fn)


def all_builtins():
    pat = re.compile(r'symtab_add_builtin\(\s*"([^"]+)"')
    out = {}
    for path in _walk_c():
        rel = os.path.relpath(path, ROOT)
        with open(path, encoding="utf-8", errors="replace") as fh:
            for m in pat.finditer(fh.read()):
                out.setdefault(m.group(1), rel)
    return out


def live_names():
    """Names["*"] from the running system.

    The C scan alone is not the inventory: `Ordering`, `Standardize` and
    everything else defined in src/internal/*.m has no symtab_add_builtin call
    to find, and a head defined in Mathilda's own language can sit off the
    buffer exactly as a C one can -- `Ordering` materialised 200000 elements
    and the scan could not even name it.  Asking the binary is the only
    complete answer.

    Option symbols, colours and pattern heads come back too.  They cost one
    discovery probe each and are then dropped, because none of them accepts an
    array shape -- which is cheaper than maintaining a list of what to exclude,
    and cannot go stale."""
    try:
        fd, path = tempfile.mkstemp(suffix=".m")
        with os.fdopen(fd, "w") as fh:
            fh.write('Scan[Print["#", #] &, Names["*"]];\n')
        p = subprocess.run([MATHILDA, "-file", path], capture_output=True,
                           text=True, timeout=60)
        os.unlink(path)
        return {ln[1:].strip() for ln in p.stdout.splitlines()
                if ln.startswith("#") and ln[1:].strip()}
    except Exception:
        return set()


# ---------------------------------------------------------------------------
# Heads this sweep must not call.
#
# Two kinds, and the distinction matters.  SKIP_UNSAFE would damage the machine
# or the session -- the sweep cannot ask the question at all.  SKIP_EXPLOSIVE
# answers correctly but with an output super-linear in the input, so timing it
# at n = 200000 measures the size of the answer rather than the speed of the
# path; these are timed at the small-vector shape instead, never the large one.
# Neither list is an exemption from HAVING a fast path -- it is only a statement
# about what this particular probe can measure.
# ---------------------------------------------------------------------------

SKIP_UNSAFE = {
    "Quit", "Exit", "Abort", "Run", "RunProcess", "Pause",
    "OpenRead", "OpenWrite", "OpenAppend", "Close", "Read", "ReadLine",
    "ReadList", "ReadString", "Write", "WriteString", "WriteLine",
    "Get", "Put", "PutAppend", "Save", "Import", "Export",
    "DeleteFile", "CopyFile", "RenameFile", "CreateDirectory",
    "DeleteDirectory", "SetDirectory", "ResetDirectory",
    "Input", "InputString", "Interrupt", "Dialog",
    "Clear", "ClearAll", "Remove", "Unset", "Protect", "Unprotect",
    "SetAttributes", "ClearAttributes", "Begin", "End", "BeginPackage",
    "EndPackage", "Needs", "Install", "Uninstall", "LibraryFunctionLoad",
    "SeedRandom", "$Failed", "Throw", "Return", "Break", "Continue",
    "TimeConstrained", "MemoryConstrained", "Share", "ShareMemory",
    # Assignment and scoping heads hold their first argument; feeding them a
    # bound symbol rebinds the sweep's own preamble data underneath it.
    "Set", "SetDelayed", "UpSet", "UpSetDelayed", "TagSet", "TagSetDelayed",
    "AddTo", "SubtractFrom", "TimesBy", "DivideBy", "Increment", "Decrement",
    "PreIncrement", "PreDecrement", "AppendTo", "PrependTo", "AssociateTo",
    "Module", "Block", "With", "Function",
    # Windowing / rendering: opens a window under USE_GRAPHICS.
    "Show", "Graphics", "Plot", "ListPlot", "ListLinePlot", "ParametricPlot",
    "Plot3D", "ContourPlot", "DensityPlot", "PolarPlot", "LogPlot",
    "LogLogPlot", "LogLinearPlot", "Histogram", "BarChart", "PieChart",
    "Spectrogram", "Periodogram",
}


def _is_rendering(name):
    """Graphics heads by shape, so a new one does not have to be enumerated.

    ParametricPlot3D was missing from the list above on the first full run;
    matching the suffix is what stops the next one from being missing too."""
    return (name.endswith(("Plot", "Plot3D", "Chart", "Graphics"))
            or name.startswith(("Plot", "ListPlot", "Graphics")))

SKIP_EXPLOSIVE = {
    # Output quadratic or worse in the input length.
    "Subsets", "Tuples", "Permutations", "Outer", "Distribute",
    "KroneckerProduct", "Nearest", "Minors", "Adjugate",
    "IntegerPartitions", "Compositions", "Multinomial",
    "DiagonalMatrix", "HankelMatrix", "ToeplitzMatrix", "VandermondeMatrix",
    "HilbertMatrix", "FourierMatrix", "FourierDCTMatrix", "FourierDSTMatrix",
    "IdentityMatrix", "UnitVector", "Tensor", "TensorProduct",
    # Exact-integer growth: the ELEMENTS explode, not their count.  Prime[k]
    # and k! over Range[200000] measure GMP, not the buffer.
    "Factorial", "Factorial2", "Fibonacci", "LucasL", "Power", "Prime",
    "PartitionsP", "PartitionsQ", "BernoulliB", "EulerE", "CatalanNumber",
    "Hyperfactorial", "BarnesG", "Subfactorial", "PrimePi", "NextPrime",
    "RandomPrime", "FactorInteger", "Divisors", "DivisorSum", "Pochhammer",
    # FactorialPower[k, m] = k(k-1)...(k-m+1) is m terms; over Range[n] with m up
    # to ~n it is the same GMP blow-up as Pochhammer, so it is priced by the
    # integers, not the buffer.
    "FactorialPower",
    # Symbolic engines: correct on a list, and priced by the algebra.
    "Solve", "Reduce", "FindInstance", "Eliminate", "GroebnerBasis",
    "Integrate", "DSolve", "RSolve", "FullSimplify", "Simplify",
    "Factor", "FactorList", "Series", "Limit", "Residue", "Resultant",
    "Expand", "ExpandAll", "Together", "Apart", "Cancel", "Collect",
    "PolynomialGCD", "PolynomialLCM", "Roots", "NRoots", "NSolve",
    "RootReduce", "ToRadicals", "Discriminant", "Interpolation",
    "ListInterpolation", "InterpolatingPolynomial", "Fit", "FindFit",
    "LatticeReduce", "FindIntegerNullVector", "Orthogonalize",
    # More symbolic engines, surfaced by the first --gate-only run at the n the
    # makefile actually uses (5000) rather than the 50000 the baseline below was
    # first recorded at.  Every one is degenerate or a no-op on a machine array
    # -- D[v, 2] is 0, Variables[v] is {}, Coefficient[v, w] is a zero vector,
    # TrigExpand/PowerExpand/ComplexExpand thread as the identity, PolynomialQ
    # declines -- so the probe measures the algebra (or the Listable thread of a
    # no-op), never a fast path a buffer could take.  They belong here beside
    # Simplify and Factor, not in OFF_BUFFER, which claims a fast path is owed.
    "D", "Dt", "Coefficient", "CoefficientList", "Decompose",
    "FactorTermsList", "MinimalPolynomial", "PolynomialSqrt", "PolynomialMod",
    "PolynomialQ", "IrreduciblePolynomialQ", "Variables", "PossibleZeroQ",
    "SimplifyCount", "TrigExpand", "TrigFactor", "TrigReduce", "TrigToExp",
    "ComplexExpand", "PowerExpand",
    # Iterate-to-convergence: a doubling step over a float vector converges
    # only via overflow to Infinity, which is thousands of whole-array passes.
    "FixedPoint", "FixedPointList", "NestWhile", "NestWhileList",
}

# Heads that read an argument as a DIMENSION SPEC.  Handing them a 200000-long
# vector asks for an array of prod(1..200000) elements, which is not a slow
# probe but a fatal one -- there is no timeout short enough to make it safe, so
# these are excluded rather than shrunk.  None of them is a fast-path question:
# each is a producer, and whether its OUTPUT packs is what nd_surface_audit.py
# --survival already asks.
SKIP_DIMS = {
    "ConstantArray", "Array", "ArrayReshape", "Table", "SparseArray",
    "RandomReal", "RandomInteger", "RandomComplex", "RandomVariate",
    "RandomChoice", "RandomSample", "Range", "Subdivide", "ArrayPad",
    # PadLeft[list, n] takes n as a full dimension spec, so PadLeft[vi, wi]
    # over two 200000-long integer vectors asks for an array of prod(wi)
    # elements.  This is the one that hung the first full run.
    "PadLeft", "PadRight",
}


# ---------------------------------------------------------------------------
# Call shapes.  Discovery tries every one of these at n = 8 and keeps the ones
# that evaluate; measurement then runs the survivors at full size.
#
# The order is the priority order: when a head accepts several shapes the
# report leads with the earliest, which is the plainest "here is an array,
# do your thing" spelling.
# ---------------------------------------------------------------------------

SHAPES = [
    ("v",     "%s[vf]",            "vec"),      # H[floatvector]
    ("vi",    "%s[vi]",            "vec"),      # H[intvector]
    ("vv",    "%s[vf, wf]",        "vec"),      # H[floatvec, floatvec]
    ("vvi",   "%s[vi, wi]",        "vec"),      # H[intvec, intvec]
    ("vk",    "%s[vf, 2]",         "vec"),      # H[floatvec, int]
    ("vki",   "%s[vi, 2]",         "vec"),      # H[intvec, int]
    ("kv",    "%s[2, vf]",         "vec"),      # H[int, floatvec]
    ("m",     "%s[mf]",            "mat"),      # H[floatmatrix]
    ("mi",    "%s[mi]",            "mat"),      # H[intmatrix]
    ("fv",    "%s[2. # &, vf]",    "vec"),      # H[function, floatvec]
    ("pv",    "%s[# > 0.5 &, vf]", "vec"),      # H[predicate, floatvec]
    ("vp",    "%s[vf, # > 0.5 &]", "vec"),      # H[floatvec, predicate]
]

SHAPE_BY_ID = {s[0]: s for s in SHAPES}


def preamble(n, mdim):
    """Deterministic data, identical under every environment the sweep uses."""
    return f"""
n = {n};
vf = Range[1., n]/n;
wf = Reverse[vf];
vi = Range[n];
wi = Mod[Range[n]*7919, 1000] + 1;
mf = Table[If[i == j, {mdim}., 1./(1. + Abs[i - j])], {{i, {mdim}}}, {{j, {mdim}}}];
mi = Table[Mod[i*j, 97] + 1, {{i, {mdim}}}, {{j, {mdim}}}];
"""


# ---------------------------------------------------------------------------
# Phase 1 -- discovery.  Which shapes does this head accept?
#
# "Accept" is `Head[result] =!= H`.  A builtin that cannot handle its argument
# returns the call unevaluated (SPEC.md §13: NULL means "I can't evaluate
# this"), so the head coming back unchanged is the system's own report that the
# shape was rejected.  The few heads that legitimately answer with their own
# head are listed in SELF_HEADED and judged on the result's LENGTH instead.
# ---------------------------------------------------------------------------

# `Head[result] =!= H` alone is FAR too permissive, and the first gate run was
# useless because of it: almost every head is Listable, so `BesselJ[v]` threads
# into a List of 50000 UNEVALUATED BesselJ calls -- whose head is List, so the
# probe read it as accepted.  Half the symbol table came back as a finding.
#
# The test is therefore "the head does not appear anywhere in the answer"
# (FreeQ), which rejects a threaded list of unevaluated calls and keeps
# Sort[v], Total[v], Chop[v] and every genuine answer.  It also correctly
# rejects Sin[intVector], which is symbolic in the interpreter and SHOULD
# degrade rather than compute in float.
SELF_HEADED = {"List", "Complex", "Interval", "SparseArray", "NDArray",
               "Association", "Rule", "RuleDelayed", "Span", "Range"}

DISCOVER_TMPL = """
r = %(expr)s;
Print["@", "%(id)s", "\\t",
      If[Head[r] === %(head)s || !FreeQ[r, %(head)s], "no", "yes"], "\\t",
      ToString[Head[r]], "\\t", If[ListQ[r] || NDArrayQ[r], Length[r], -1]];
Clear[r];
"""


def discover_source(items, n, mdim):
    out = [preamble(n, mdim)]
    for head, sid, expr in items:
        out.append(DISCOVER_TMPL % {"expr": expr, "id": f"{head}/{sid}",
                                    "head": head})
    return "\n".join(out)


def parse_discover(text):
    out = {}
    for line in text.splitlines():
        if not line.startswith("@"):
            continue
        parts = line[1:].split("\t")
        if len(parts) < 4:
            continue
        out[parts[0].strip()] = {
            "ok": parts[1].strip() == "yes",
            "head": parts[2].strip(),
            "len": parts[3].strip(),
        }
    return out


# ---------------------------------------------------------------------------
# Phase 2 -- measurement.
#
# Minimum of REPS, with the maximum kept as a caching tripwire: a head that
# memoises answers the second call in ~0 s, and a minimum over a memoised
# repeat is a measurement of the cache, not of the path (see
# docs/experiments/README.md, and the `eval memo` note in tasks/lessons.md).
# A row whose max/min exceeds MEMO_RATIO is reported as suspect rather than
# silently believed.
# ---------------------------------------------------------------------------

MEMO_RATIO = 8.0

TIME_TMPL = """
tmin = 1000000.; tmax = 0.;
Do[dt = First[AbsoluteTiming[r = %(expr)s]];
   tmin = Min[tmin, dt]; tmax = Max[tmax, dt], {%(reps)d}];
Print["@", "%(id)s", "\\t", tmin, "\\t", tmax, "\\t",
      If[NDArrayQ[r], "packed", "plain"]];
Clear[r];
"""


def time_source(items, n, mdim):
    out = [preamble(n, mdim)]
    for head, sid, expr in items:
        out.append(TIME_TMPL % {"expr": expr, "id": f"{head}/{sid}",
                                "reps": REPS})
    return "\n".join(out)


def parse_time(text):
    out = {}
    for line in text.splitlines():
        if not line.startswith("@"):
            continue
        parts = line[1:].split("\t")
        if len(parts) < 4:
            continue
        try:
            tmin = float(parts[1])
            tmax = float(parts[2])
        except ValueError:
            continue
        out[parts[0].strip()] = {"min": tmin, "max": tmax,
                                 "packed": parts[3].strip() == "packed"}
    return out


# ---------------------------------------------------------------------------
# Runner.  One process per batch, and a batch that blows the timeout is retried
# one probe per process so a single hang costs one row rather than the batch.
# ---------------------------------------------------------------------------

def run_batch(items, source_fn, parse_fn, n, mdim, env_extra=None, tmpdir=None,
              timeout=None):
    if not items:
        return {}
    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)
    fd, path = tempfile.mkstemp(suffix=".m", dir=tmpdir)
    with os.fdopen(fd, "w") as fh:
        fh.write(source_fn(items, n, mdim))
    try:
        proc = subprocess.run([MATHILDA, "-file", path], env=env,
                              capture_output=True, text=True,
                              timeout=timeout or TIMEOUT)
        return parse_fn(proc.stdout)
    except subprocess.TimeoutExpired:
        return None
    finally:
        os.unlink(path)


def run_all(items, source_fn, parse_fn, n, mdim, batch, env_extra=None,
            label="", timeout=None):
    got = {}
    with tempfile.TemporaryDirectory() as td:
        i = 0
        while i < len(items):
            chunk = items[i:i + batch]
            res = run_batch(chunk, source_fn, parse_fn, n, mdim, env_extra, td,
                            timeout)
            if res is None and len(chunk) > 1:
                # BISECT.  Retrying one probe at a time is what made a large
                # batch expensive: one bad probe in 128 cost 128 re-runs, each
                # paying the timeout's worth of nothing.  Halving costs
                # ~2*log2(n) runs and isolates the same probe.
                stack = [chunk]
                while stack:
                    part = stack.pop()
                    r1 = run_batch(part, source_fn, parse_fn, n, mdim,
                                   env_extra, td, timeout)
                    if r1:
                        got.update(r1)
                    elif len(part) > 1:
                        mid = len(part) // 2
                        stack.append(part[:mid])
                        stack.append(part[mid:])
                    # len(part) == 1 and no result: that probe does not answer.
            elif res:
                got.update(res)
            i += batch
            if label:
                sys.stderr.write(f"\r  {label}: {min(i, len(items))}/"
                                 f"{len(items)}   ")
                sys.stderr.flush()
    if label:
        sys.stderr.write("\n")
    return got


# ---------------------------------------------------------------------------
# The gate diagnostic: which heads made the transparency gate materialise?
# Confirmation for a row the clock has already flagged, and the only one of the
# three signals that names the mechanism.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# The GATE PASS -- the primary detector.
#
# MATHILDA_PACK_DIAG=gate makes pack.c count, per head, how many elements the
# transparency gate materialised: one boxed Expr per element, built before the
# builtin ever ran.  That is the exact mechanism a missing fast path shows up
# as, and it is a COUNT rather than a duration -- deterministic, immune to what
# else is running on the machine, and complete in minutes where timing every
# head twice takes hours.
#
# The timing half is still worth having, but it is the SEVERITY, not the
# detector: run it over the heads this flags, on an idle machine.
#
# Two things have to be read together, which is why one pass collects both:
#   * the gate count says the buffer was thrown away;
#   * `evaluated` says the call meant anything.  A head that returns its
#     argument unevaluated ALSO trips the post-gate as the node comes to rest,
#     and that is not a finding -- it is a head being handed a shape it does
#     not take.
# ---------------------------------------------------------------------------

GATE_TMPL = """
r = %(expr)s;
Print["@", "%(id)s", "\\t", If[Head[r] === %(head)s, "no", "yes"]];
Clear[r];
"""


def gate_source(items, n, mdim):
    out = [preamble(n, mdim)]
    for head, sid, expr in items:
        out.append(GATE_TMPL % {"expr": expr, "id": f"{head}/{sid}",
                                "head": head})
    return "\n".join(out)


GATE_ROW = re.compile(r"^\s{7}(\S+)\s+(\d+)\s+(\d+)\s*$")


def gate_pass(items, n, mdim, batch=48, timeout=60):
    """head -> {'elements': int, 'calls': int}, plus which ids evaluated.

    Batches probes and, on a batch timeout, keeps whatever the batch printed
    before it (partial). This is FAST -- a release-gate run stays in the tens of
    minutes -- but it means a batch that times out loses every probe after its
    first slow one, so the counted set flaps by a handful of near-timeout heads
    from run to run. The ratchet tolerates this by keying its exit status on NEW
    heads only and by tracking a deliberately GENEROUS OFF_BUFFER superset, so a
    flap shows up at worst as an advisory "FIXED" line, never a failure. A
    bisecting version that drops only genuinely-hung probes was tried and
    reverted: with `timeout` at 60 s it re-ran the slow batches down to single
    probes and pushed one run past two hours, which is the wrong trade for a
    gate whose whole value is being minutes rather than the timing pass's
    hours."""
    totals, evaluated = {}, {}
    env = dict(os.environ)
    env["MATHILDA_PACK_DIAG"] = "gate"
    with tempfile.TemporaryDirectory() as td:
        for i in range(0, len(items), batch):
            chunk = items[i:i + batch]
            fd, path = tempfile.mkstemp(suffix=".m", dir=td)
            with os.fdopen(fd, "w") as fh:
                fh.write(gate_source(chunk, n, mdim))
            try:
                p = subprocess.run([MATHILDA, "-file", path], env=env,
                                   capture_output=True, text=True,
                                   timeout=timeout)
                out, err = p.stdout, p.stderr
            except subprocess.TimeoutExpired as e:
                out = e.stdout.decode() if e.stdout else ""
                err = e.stderr.decode() if e.stderr else ""
            finally:
                os.unlink(path)
            for line in out.splitlines():
                if line.startswith("@"):
                    parts = line[1:].split("\t")
                    if len(parts) >= 2:
                        evaluated[parts[0].strip()] = parts[1].strip() == "yes"
            for line in err.splitlines():
                m = GATE_ROW.match(line)
                if m:
                    t = totals.setdefault(m.group(1),
                                          {"calls": 0, "elements": 0})
                    t["calls"] += int(m.group(2))
                    t["elements"] += int(m.group(3))
            sys.stderr.write(f"\r  gate: {min(i + batch, len(items))}/"
                             f"{len(items)}   ")
            sys.stderr.flush()
    sys.stderr.write("\n")
    return totals, evaluated


# ---------------------------------------------------------------------------
# The RATCHET.
#
# These heads accept an array shape and materialise it: real work, not yet
# done. Listing them is what lets the exit status mean "something NEW fell off
# the buffer" instead of "the backlog is still there". Every entry is a
# candidate fast path; several are cheap.
#
# Recorded with the EXACT command the makefile runs -- `nd_fastpath_sweep.py
# --gate-only`, i.e. n=5000.  The first baseline (2026-08-02) was taken at
# --n 50000, where dozens of heads time out in the gate pass and never get
# classified, so it silently under-recorded the backlog; keying the record to a
# different n than the check runs at is what let this list go stale unnoticed.
# Last re-recorded 2026-08-03.
# ---------------------------------------------------------------------------

OFF_BUFFER = {
    # Element-wise, and cheap: the answer is one value per element and the input
    # is already a buffer.  Chop is the clearest -- a threshold compare.
    "Chop", "Im", "Rationalize", "Numerator", "Denominator",
    "MantissaExponent", "RealDigits",
    # Integer predicates and number-theoretic maps over an int64 buffer.  Several
    # already have ndinteger.c neighbours (EulerPhi, MoebiusMu) to copy.
    "EvenQ", "OddQ", "CoprimeQ", "PrimeNu", "PrimeOmega", "LiouvilleLambda",
    "JacobiSymbol", "DigitSum", "IntegerExponent",
    # Search and selection.  Position and Count are the two a numeric workload
    # reaches most; both are a single pass over the buffer.
    "Position", "Count", "MemberQ", "Cases", "FirstCase", "Delete",
    "ReplacePart", "Level", "Apply", "Thread",
    # Grouping and ordering -- a sort or a hash over machine words, the same
    # shape as the Tally/DeleteDuplicates work already done in ndreduce.c.
    "Split", "Gather", "GatherBy", "PositionIndex", "Lookup", "OrderedQ",
    "MinimalBy", "ReverseSortBy",
    # Shape predicates: they read the dims, or at most one pass, and the buffer
    # already knows its shape.
    "SquareMatrixQ", "SymmetricMatrixQ", "DiagonalMatrixQ",
    "UpperTriangularMatrixQ", "ArrayFlatten",
    # Comparisons and same-tests.  The answer is a Boolean or a boolean List, so
    # this buys the INPUT side only -- there is no bool dtype (§13 gap C.1).
    "Less", "LessEqual", "Greater", "GreaterEqual", "Inequality",
    "Unequal", "SameQ", "UnsameQ",
    # Two ARRAY operands.  ndarray_map_binary broadcasts one array against a
    # scalar and has no array-against-array map, so these need that first.
    # BesselY joined its J/I/K siblings from the n=5000 run.
    "BesselJ", "BesselI", "BesselK", "BesselY", "QuotientRemainder",
    # ------------------------------------------------------------------
    # Everything below was recorded 2026-08-03 from the n=5000 gate pass -- the
    # backlog the n=50000 baseline could not see.  None touches runtime; each is
    # a head that materialises a packed argument it accepts and could later take
    # a buffer path (or, for the pass-throughs, join pack.c's AWARE list).
    #
    # Elementwise numeric / number-theoretic maps -- one value per element over a
    # buffer, the same shape as the Chop/PrimeNu group above.
    "PrimeQ", "SquareFreeQ", "Divisible", "ReIm", "RealExponent",
    "PrimitiveRoot", "PrimitiveRootList", "ContinuedFraction",
    # Precision control.  Precision/Accuracy were pulled from AWARE (an exact
    # element has no machine precision), so they materialise by design; SetPrecision
    # and SetAccuracy rebuild every element at a new precision.
    "Precision", "Accuracy", "SetPrecision", "SetAccuracy",
    # Type / identity predicates.  The answer does not depend on the elements at
    # all -- a packed List is not an Integer, a Number or True -- so these
    # materialise only because the gate fires before them; each is really an
    # AWARE pass-through waiting to be proven safe.
    "IntegerQ", "NumberQ", "NumericQ", "MachineNumberQ", "TrueQ",
    "AssociationQ", "DirectedGraphQ", "GraphQ", "ListQ",
    # Selection and grouping keyed by a function or a pattern -- siblings of the
    # MinimalBy / Cases / DeleteDuplicates entries already on the buffer.
    "MaximalBy", "DeleteCases", "DeleteDuplicatesBy", "CountsBy", "GroupBy",
    "Pick", "DeleteMissing",
    # The rule / pattern engine.  Replace and friends walk the materialised tree;
    # a packed numeric argument is boxed whole before a single rule is tried.
    "Replace", "ReplaceAll", "ReplaceList", "ReplaceRepeated", "MatchQ",
    # Association construction and key queries over a machine vector.
    "AssociationMap", "AssociationThread", "KeyExistsQ", "KeyFreeQ",
    "KeyMemberQ",
    # Canonical comparison of two whole arrays -- could compare buffers directly
    # rather than materialising both into Expr trees.
    "Order",
    # Functional / meta heads that only move the value around or format it, so a
    # packed argument survives untouched -- pass-through candidates for AWARE.
    # Identity is the purest of these; Goto passes its target through untouched.
    "Composition", "ComposeList", "Through", "ReleaseHold", "Sequence",
    "Options", "Assuming", "ToString", "IntegerString", "Identity", "Goto",
    # String heads handed a numeric vector: they format or split it, boxing every
    # element to a digit string first.  StringSplit/StringTrim are degenerate on
    # numbers but still materialise; StringRiffle joins the whole vector.
    "StringRiffle", "StringSplit", "StringTrim",
    # Graph constructor over a coordinate/edge vector.
    "PathGraph",
}


def gate_main(args, cand):
    n = args.n
    mdim = max(4, int(n ** 0.5))

    # Phase 1: discovery at n = 8, which is nearly hang-free and cheap.
    all_probes = [(h, sid, tmpl % h) for h in cand for sid, tmpl, _ in SHAPES
                  if h not in SKIP_EXPLOSIVE]
    sys.stderr.write(f"discovering {len(all_probes)} shapes over {len(cand)} "
                     f"heads\n")
    disc = run_all(all_probes, discover_source, parse_discover, 8, 4,
                   DISCOVER_BATCH, label="discover", timeout=DISCOVER_TIMEOUT)

    # Phase 2: the gate diagnostic over ONLY the shapes that mean something.
    # Running every shape and reading a per-head TOTAL mixes in the calls that
    # came back unevaluated -- the gate fires for those too -- so the count
    # would not belong to the call the report names.
    probes = [(h, sid, e) for h, sid, e in all_probes
              if (disc.get(f"{h}/{sid}") or {}).get("ok")]
    shapes_ok = {}
    for h, sid, expr in probes:
        shapes_ok.setdefault(h, []).append(expr)
    sys.stderr.write(f"gate pass: {len(probes)} accepted shapes over "
                     f"{len(shapes_ok)} heads at n={n}\n")
    totals, _ = gate_pass(probes, n, mdim)

    # Phase 3: VERIFY.  The pass above runs 16 probes per process and the gate
    # counter is per head over the whole process, so a head can be credited
    # with a materialisation that belongs to a neighbour -- the first full run
    # named Exp, Gamma, Variance, Clip and thirty more that do not materialise
    # at all when asked on their own.  Re-ask the survivors in ONE process,
    # one explicit expression per head, and keep only what answers again.
    #
    # The size drops here on purpose: "did the gate fire" is a BOOLEAN, so a
    # small vector answers it, and at 50000 the quadratic shapes (Position and
    # Count over two vectors) never finish.
    # A head the diagnostic names need not be one this run PROBED: an internal
    # `List[...]` assembling a result materialises too, and the counter records
    # it by name like any other. Those have no probe expression to re-ask with,
    # and they are not findings about a builtin -- drop them before the verify
    # pass rather than crash on the missing key, which is what the first
    # baseline-validation run did (KeyError: 'List').
    totals = {h: t for h, t in totals.items() if h in shapes_ok}

    if totals:
        cand = sorted(totals, key=lambda h: -totals[h]["elements"])
        vn, vm = 3000, 54
        lines = [preamble(vn, vm)]
        for h in cand:
            lines.append(f"q = {shapes_ok[h][0]};\nClear[q];")
        fd, path = tempfile.mkstemp(suffix=".m")
        with os.fdopen(fd, "w") as fh:
            fh.write("\n".join(lines) + '\nPrint["done"];\n')
        env = dict(os.environ)
        env["MATHILDA_PACK_DIAG"] = "gate"
        confirmed = {}
        try:
            pr = subprocess.run([MATHILDA, "-file", path], env=env,
                                capture_output=True, text=True, timeout=600)
            for line in pr.stderr.splitlines():
                mm = GATE_ROW.match(line)
                if mm and mm.group(1) in totals:
                    confirmed[mm.group(1)] = {"calls": int(mm.group(2)),
                                              "elements": int(mm.group(3))}
        except subprocess.TimeoutExpired:
            sys.stderr.write("verify pass timed out; reporting unverified\n")
            confirmed = totals
        finally:
            os.unlink(path)
        sys.stderr.write(f"verify: {len(confirmed)}/{len(totals)} heads "
                         f"materialise when asked on their own\n")
        totals = confirmed

    rows = []
    for h, t in totals.items():
        if h not in shapes_ok:
            continue                     # nothing it accepts; not a finding
        rows.append(dict(head=h, calls=t["calls"], elements=t["elements"],
                         shapes=shapes_ok[h]))
    rows.sort(key=lambda r: -r["elements"])

    print(f"\n{'head':<28} {'calls':>7} {'elements':>14}  a shape it accepts")
    print("-" * 92)
    for r in rows:
        print(f"{r['head']:<28} {r['calls']:>7} {r['elements']:>14}  "
              f"{r['shapes'][0]}")
    print(f"\n{len(rows)} heads made the transparency gate materialise a packed "
          f"argument they accept.\nEach is a fast path not taken: the buffer was "
          f"turned into one Expr per element\nbefore the builtin ran. Time them "
          f"with --only <head> to rank the severity.")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(rows, fh, indent=1)
        print(f"wrote {args.json}")

    # The ratchet.  Most of the standing list is real work that has not been
    # done yet, and a gate that fails on all of it from the day it lands is
    # noise -- nobody can tell a new regression from the backlog.  So the exit
    # status keys on the DIFFERENCE: a head that materialises and is not in
    # OFF_BUFFER fails, and a head in OFF_BUFFER that has stopped is reported
    # so the line gets deleted.  A --only run never fails on absences.
    names = {r["head"] for r in rows}
    new_heads = sorted(names - OFF_BUFFER)
    if new_heads:
        print("\nNEW: these heads materialise a packed argument they accept, "
              "and are not in OFF_BUFFER:")
        for h in new_heads:
            print(f"  {h}")
        print("Give the head a buffer path and put it on pack.c's AWARE list, "
              "or add it to\nOFF_BUFFER in "
              f"{os.path.relpath(__file__, ROOT)} to record it as known.")
        return 1
    if not args.only:
        fixed = sorted(OFF_BUFFER - names)
        if fixed:
            print("\nFIXED since the baseline was recorded — delete these from "
                  "OFF_BUFFER:")
            for h in fixed:
                print(f"  {h}")
        print(f"\nOK: no new heads off the buffer ({len(names & OFF_BUFFER)} "
              f"known, tracked in OFF_BUFFER).")
    return 0


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=None,
                    help="elements in the probe vectors "
                         "(default 100000 timing, 5000 gate)")
    ap.add_argument("--only", default="", help="comma-separated heads")
    ap.add_argument("--discover-only", action="store_true")
    ap.add_argument("--all", action="store_true", help="print every measured row")
    ap.add_argument("--json", default="")
    ap.add_argument("--ns-threshold", type=float, default=40.0,
                    help="ns/elem above which a row is a candidate (default 40)")
    ap.add_argument("--ratio-threshold", type=float, default=1.6,
                    help="t_nopack/t_packed below which packing bought nothing")
    ap.add_argument("--gate", action="store_true",
                    help="confirm flagged rows with MATHILDA_PACK_DIAG=gate")
    ap.add_argument("--gate-only", action="store_true",
                    help="just the gate pass: which heads materialise a buffer "
                         "they accept? deterministic and minutes, not hours")
    args = ap.parse_args()

    if not os.path.exists(MATHILDA):
        sys.exit(f"no binary at {MATHILDA}; run make first")

    heads = all_builtins()
    # Union with the live symbol table: the C scan cannot see a head defined in
    # src/internal/*.m, and those can be off the buffer too.
    for nm in live_names():
        heads.setdefault(nm, "(internal .m)")
    wanted = [h.strip() for h in args.only.split(",") if h.strip()]
    if wanted:
        cand = [h for h in wanted if h in heads] or wanted
    else:
        cand = sorted(h for h in heads
                      if h not in SKIP_UNSAFE and h not in SKIP_DIMS
                      and not _is_rendering(h) and not h.startswith("$")
                      and "`" not in h)

    if args.gate_only:
        # "Did the gate fire" is a BOOLEAN, so the gate phase does not need a
        # big vector -- it needs one above the packing threshold (4).  The
        # first full run used 50000 and spent its time building a 223x223
        # preamble matrix 118 times over for an answer that does not depend on
        # the size.  Timing is the mode that needs the elements.
        if args.n is None:
            args.n = 5000
        return gate_main(args, cand)

    if args.n is None:
        args.n = 100000
    mdim = max(4, int(args.n ** 0.5))
    small_n, small_mdim = 8, 4

    # --- phase 1: which shapes does each head accept? ----------------------
    probes = [(h, sid, tmpl % h) for h in cand for sid, tmpl, _ in SHAPES]
    sys.stderr.write(f"discovering {len(probes)} shapes over {len(cand)} heads\n")
    disc = run_all(probes, discover_source, parse_discover, small_n, small_mdim,
                   DISCOVER_BATCH, label="discover", timeout=DISCOVER_TIMEOUT)

    accepted = {}
    for h, sid, expr in probes:
        d = disc.get(f"{h}/{sid}")
        if d and d["ok"]:
            accepted.setdefault(h, []).append((sid, expr, d))

    if args.discover_only:
        for h in sorted(accepted):
            shapes = ", ".join(s for s, _, _ in accepted[h])
            print(f"{h:<28} {shapes}")
        print(f"\n{len(accepted)}/{len(cand)} heads accept at least one array shape")
        return 0

    # Time at most two shapes per head: the plainest accepted one, and the
    # plainest int64 one if it differs.  More than that measures the same path
    # repeatedly and lengthens the sweep for nothing.
    to_time = []
    for h in sorted(accepted):
        picks, seen_int = [], False
        for sid, expr, d in accepted[h]:
            is_int = sid in ("vi", "vvi", "vki", "mi")
            if not picks:
                picks.append((sid, expr))
                seen_int = is_int
            elif is_int and not seen_int:
                picks.append((sid, expr))
                break
        n_here = args.n if h not in SKIP_EXPLOSIVE else SMALL_N
        for sid, expr in picks:
            to_time.append((h, sid, expr, n_here))

    # --- phase 2: measure, packed and unpacked -----------------------------
    big = [(h, sid, e) for h, sid, e, nn in to_time if nn == args.n]
    small = [(h, sid, e) for h, sid, e, nn in to_time if nn != args.n]

    sys.stderr.write(f"timing {len(to_time)} probes at n={args.n}\n")
    tp = run_all(big, time_source, parse_time, args.n, mdim, TIME_BATCH,
                 label="packed")
    tu = run_all(big, time_source, parse_time, args.n, mdim, TIME_BATCH,
                 {"MATHILDA_NO_PACK": "1"}, label="nopack")
    if small:
        tp.update(run_all(small, time_source, parse_time, SMALL_N, SMALL_MDIM,
                          TIME_BATCH, label="packed(sm)"))
        tu.update(run_all(small, time_source, parse_time, SMALL_N, SMALL_MDIM,
                          TIME_BATCH, {"MATHILDA_NO_PACK": "1"},
                          label="nopack(sm)"))

    rows = []
    for h, sid, expr, nn in to_time:
        key = f"{h}/{sid}"
        p, u = tp.get(key), tu.get(key)
        if not p:
            rows.append(dict(head=h, shape=sid, expr=expr, n=nn,
                             status="no-complete"))
            continue
        kind = SHAPE_BY_ID[sid][2]
        elems = nn if kind == "vec" else (mdim * mdim if nn == args.n
                                          else SMALL_MDIM * SMALL_MDIM)
        ns = p["min"] / elems * 1e9
        ratio = (u["min"] / p["min"]) if (u and p["min"] > 0) else float("nan")
        memo = (p["max"] / p["min"]) if p["min"] > 0 else 1.0
        rows.append(dict(head=h, shape=sid, expr=expr, n=nn, elems=elems,
                         ns_per_elem=ns, t_packed=p["min"], t_nopack=
                         (u["min"] if u else None), ratio=ratio,
                         packed_result=p["packed"], memo=memo,
                         status="ok"))

    flagged = [r for r in rows
               if r["status"] == "ok"
               and r["ns_per_elem"] >= args.ns_threshold
               and (r["ratio"] != r["ratio"] or r["ratio"] < args.ratio_threshold)]
    flagged.sort(key=lambda r: -r["ns_per_elem"])

    if args.gate and flagged:
        sys.stderr.write("confirming flagged rows under the gate diagnostic\n")
        # One probe per process, so the per-head total belongs to this row.
        for r in flagged:
            t, _ = gate_pass([(r["head"], r["shape"], r["expr"])], args.n, mdim,
                             batch=1)
            ent = t.get(r["head"])
            r["gate_elems"] = ent["elements"] if ent else None

    # --- report -------------------------------------------------------------
    show = rows if args.all else flagged
    if not args.all:
        show.sort(key=lambda r: -r.get("ns_per_elem", 0))
    else:
        show = sorted(rows, key=lambda r: -r.get("ns_per_elem", 0))

    print(f"\n{'head':<26} {'shape':<5} {'ns/elem':>10} {'nopack/packed':>14} "
          f"{'result':>8}  expr")
    print("-" * 96)
    for r in show:
        if r["status"] != "ok":
            print(f"{r['head']:<26} {r['shape']:<5} {'--':>10} {'--':>14} "
                  f"{'--':>8}  {r['expr']}   ({r['status']})")
            continue
        res = "packed" if r["packed_result"] else "plain"
        note = ""
        if r["memo"] > MEMO_RATIO:
            note = "  [memo?]"
        if r.get("gate_elems"):
            note += f"  [gate {r['gate_elems']}]"
        print(f"{r['head']:<26} {r['shape']:<5} {r['ns_per_elem']:>10.1f} "
              f"{r['ratio']:>14.2f} {res:>8}  {r['expr']}{note}")

    ok = [r for r in rows if r["status"] == "ok"]
    print(f"\nmeasured {len(ok)} probes over {len(accepted)} heads; "
          f"{len(flagged)} flagged "
          f"(>= {args.ns_threshold:g} ns/elem and packing bought < "
          f"{args.ratio_threshold:g}x)")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"n": args.n, "rows": rows}, fh, indent=1)
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
