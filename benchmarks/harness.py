#!/usr/bin/env python3
"""benchmarks/harness.py -- the ONLY copy of the reporting helpers.

Imported by every experiment with two lines of preamble::

    import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from harness import bench, check, require, bench_if

Python has ``__file__``, so unlike the ``.m`` side this is cwd-independent: a
``.py`` experiment runs by hand from anywhere.  (The ``.m`` side must be run
from inside its folder because Mathilda's ``Get`` resolves against cwd and
``$InputFileName`` does not exist.)

THE EMIT CONTRACT -- one tagged tab-separated line per CASE (one named kernel,
timed in every system and joined across them by its label), identical to
``harness.m`` so ``run_all.py`` has a single parser for all three systems::

    BENCH   <label>  <milliseconds>
    CHECK   <label>  <value>
    REQUIRE <head>   present|absent
    SKIP    <label>  <head>

``<label>`` is the join key between the ``.m`` and ``.py`` halves of an
experiment.  ``run_all.py`` fails an experiment whose halves emit different
label sets -- silent case-drift is how a table starts comparing two different
programs.
"""

import importlib
import time

# ---- timing method -------------------------------------------------------

BENCH_REPS = 3


def bench(label, fn, reps=BENCH_REPS):
    """One untimed warm-up, then the MINIMUM of `reps` timed runs.

    The minimum, not the mean: we are measuring the cost of the work, and every
    source of noise on a loaded machine can only add.  ``perf_counter``, which
    is the wall clock -- matching ``AbsoluteTiming`` on the CAS side rather
    than ``Timing[]``, which sums CPU time over threads.
    """
    fn()                                          # warm-up, untimed
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    print("BENCH\t%s\t%.3f" % (label, 1000.0 * min(ts)), flush=True)


def bench_once(label, fn):
    """One timed run, no warm-up -- for cases already expensive on their own.

    A 10**6-digit pi or a 1050**2 matmul x12 gains nothing from a warm-up and
    would take twice as long with one.  Used by 31-wolframmark.
    """
    t0 = time.perf_counter()
    fn()
    print("BENCH\t%s\t%.3f" % (label, 1000.0 * (time.perf_counter() - t0)),
          flush=True)


# ---- the value gate ------------------------------------------------------


def check(label, value):
    """A cheap scalar recorded from every system and compared by run_all.py.

    A timing case is MEANINGLESS until the systems agree on the check: the check
    pins the algorithm, so the columns cannot silently be timing three
    different computations.  A case whose check disagrees is reported CHECK-FAIL
    and its timing is DISCARDED rather than shown with a caveat.
    """
    print("CHECK\t%s\t%s" % (label, value), flush=True)


# ---- the coverage / absence register -------------------------------------


def require(spec):
    """Declare a dependency and report whether this system has it.

    `spec` is ``"module"``, ``"module:attr"``, or a list of either.  Feeds the
    coverage percentage and ABSENT.md.  Kept symmetrical with the ``.m`` side's
    ``require`` so one parser handles both, even though a Python dependency is
    essentially always present -- the asymmetry IS the finding: what Mathilda
    lacks and sympy/scipy have is exactly what the register is for.
    """
    if isinstance(spec, (list, tuple)):
        for s in spec:
            require(s)
        return
    mod, _, attr = spec.partition(":")
    try:
        m = importlib.import_module(mod)
        ok = (not attr) or hasattr(m, attr)
    except Exception:
        ok = False
    print("REQUIRE\t%s\t%s" % (spec, "present" if ok else "absent"), flush=True)


def have(spec):
    mod, _, attr = spec.partition(":")
    try:
        m = importlib.import_module(mod)
        return (not attr) or hasattr(m, attr)
    except Exception:
        return False


def bench_if(label, spec, fn, reps=BENCH_REPS):
    """Run a case only if the dependency exists; otherwise emit SKIP.

    Distinguishes "absent" from "silently never measured" -- the difference
    between an honest absence register and a table with holes in it.
    """
    if have(spec):
        bench(label, fn, reps)
    else:
        print("SKIP\t%s\t%s" % (label, spec), flush=True)


def check_if(label, spec, value_fn):
    if have(spec):
        check(label, value_fn())


# ---- deterministic seeding ----------------------------------------------


def seed():
    """Seed numpy's legacy global RNG, matching each experiment's SeedRandom[1].

    The three systems cannot be made to draw the SAME stream, which is why
    value checks use exactly-reproducible integer data or a separate small
    deterministic case -- never a float sum over random data.
    """
    import numpy as np
    np.random.seed(1)
