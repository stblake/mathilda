#!/usr/bin/env python3
"""Measure per-function timings for the documentation pages.

Writes ``site/perf.json``, which ``site/generate.py`` folds into each page's
Performance section. Kept out of the generator on purpose: timing hundreds of
functions takes minutes, while regenerating the prose takes seconds, and tying
the two together would make every documentation edit slow.

Two sources of probes:

* A GENERIC probe for anything that accepts a list of machine reals. Most of the
  list, statistics and elementary-function surface is of that shape, so this
  covers a wide sweep with no per-function work. Whether a function accepts the
  shape is decided by TRIAL against the built binary, not by a hand-maintained
  list that would drift.
* A CURATED table for functions whose interesting cost is not "a list of reals"
  -- clustering methods, integration, factorisation. These carry their own
  expression and a label saying what is being measured.

Every timing is wall-clock around a single evaluation of an already-constructed
input, so list construction is not counted. Sizes climb until a case exceeds
``--budget`` seconds, which keeps the slow functions from dominating the run.

Usage:
    python3 tools/docs_perf.py                 # everything
    python3 tools/docs_perf.py FindClusters    # named functions only
    python3 tools/docs_perf.py --budget 2.0
"""

import argparse
import json
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MATHILDA = ROOT / "Mathilda"
OUT = ROOT / "site" / "perf.json"

SIZES = [1000, 10_000, 100_000]

# Functions whose real cost is not "a list of machine reals". The expression
# uses %(n)d for the size. `label` says what the row is measuring, since for
# these the function name alone does not.
CURATED = {
    "FindClusters": [
        ("Agglomerate (default)", "FindClusters[data]", None),
        ("KMeans, 5 clusters", 'FindClusters[data, 5, Method -> "KMeans"]', None),
        # Vector input is quadratic and capped at 2,000 points, and `data` holds
        # 2n reals for n points -- so this case carries its own sizes rather
        # than being probed where the call is refused.
        ("2-D points", "FindClusters[Partition[data, 2]]", [500, 1000, 2000]),
    ],
    "Nearest": [("nearest to a target", "Nearest[data, 0.5]", None)],
    "EuclideanDistance": [("two n-vectors",
                           "EuclideanDistance[data, Reverse[data]]", None)],
}

# Setup for the curated entries: what `data` is bound to.
CURATED_SETUP = "data = RandomReal[{0, 1}, %(n)d];"


def run_script(body, timeout):
    """Evaluate a script and return its stdout, or None on failure/timeout."""
    script = ROOT / ".docs_perf_tmp.m"
    script.write_text(body)
    try:
        p = subprocess.run([str(MATHILDA), "-file", str(script)],
                           capture_output=True, text=True, timeout=timeout)
        return p.stdout
    except (subprocess.TimeoutExpired, OSError):
        return None
    finally:
        script.unlink(missing_ok=True)


def time_expression(setup, expr, head, timeout):
    """Seconds for one evaluation of `expr`, or None if it did not evaluate.

    AbsoluteTiming is the kernel's own measurement, so subprocess start-up and
    parsing are outside the number.

    The evaluated check is not a nicety. A declined call returns instantly with
    its head intact, and timing that produces a fast, plausible, meaningless
    number: measuring FindClusters on 2-D points this way reported 96 ms at
    n=1,000 but 968 us at n=10,000, because the larger input exceeded the
    2,000-point cap and was refused. A timing that falls as the input grows is
    the signature of this mistake, and the only reliable fix is to confirm the
    call produced a result at all."""
    out = run_script(
        f"{setup}\n"
        f"r = 0;\n"
        f"t = AbsoluteTiming[r = {expr};][[1]];\n"
        f"Print[If[Head[r] === {head}, \"DECLINED\", t]];\n", timeout)
    if not out:
        return None
    last = out.strip().splitlines()[-1] if out.strip() else ""
    if "DECLINED" in last:
        return None
    m = re.search(r"([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)", last)
    return float(m.group(1)) if m else None


def accepts_real_list(name, timeout):
    """Does `name[list-of-reals]` evaluate to something other than itself?

    Decided by trial: an unevaluated result comes back with the head still
    attached, which is exactly how this codebase signals "I do not handle
    this"."""
    out = run_script(
        f"d = RandomReal[{{0, 1}}, 32];\nr = {name}[d];\n"
        f"Print[If[Head[r] === {name}, \"NO\", \"YES\"]];\n", timeout)
    return bool(out) and "YES" in out


def fmt(sec):
    if sec is None:
        return None
    if sec < 1e-3:
        return f"{sec * 1e6:.0f} us"
    if sec < 1.0:
        return f"{sec * 1e3:.1f} ms"
    return f"{sec:.2f} s"


def measure(name, budget, timeout):
    """Rows for one function, or None when nothing could be measured."""
    rows = []
    cases = CURATED.get(name)
    if cases is None:
        if not accepts_real_list(name, timeout):
            return None
        cases = [("list of machine reals", f"{name}[data]", None)]

    for label, expr, sizes in cases:
        for n in (sizes or SIZES):
            setup = CURATED_SETUP % {"n": 2 * n if sizes else n}
            t = time_expression(setup, expr, name, timeout)
            if t is None:
                break                      # declined or timed out: stop growing
            rows.append({"case": label, "n": n, "time": fmt(t), "seconds": round(t, 6)})
            if t > budget:
                break                      # already slow; larger tells us nothing new
    return rows or None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*", help="functions to measure (default: all)")
    ap.add_argument("--budget", type=float, default=3.0,
                    help="stop growing n once a case exceeds this many seconds")
    ap.add_argument("--timeout", type=float, default=60.0)
    args = ap.parse_args()

    if not MATHILDA.exists():
        sys.exit("build ./Mathilda first")

    if args.names:
        names = args.names
    else:
        manifest = ROOT / "site" / "docs" / "assets" / "builtins.json"
        names = [b["name"] for b in json.loads(manifest.read_text())]

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                            capture_output=True, text=True).stdout.strip()
    host = f"{platform.machine()} {platform.system()}"

    # Merge rather than overwrite, so measuring one function does not discard
    # every other function's numbers.
    data = {}
    if OUT.exists():
        try:
            data = json.loads(OUT.read_text())
        except ValueError:
            data = {}

    measured = 0
    for i, name in enumerate(names, 1):
        rows = measure(name, args.budget, args.timeout)
        if rows:
            data[name] = {"host": host, "commit": commit, "rows": rows}
            measured += 1
        if len(names) > 20 and i % 25 == 0:
            print(f"  {i}/{len(names)} probed, {measured} measured", flush=True)

    OUT.write_text(json.dumps(data, indent=1, sort_keys=True))
    print(f"{measured} function(s) measured -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
