#!/usr/bin/env python3
"""FeynCalc microbenchmarks: Mathilda against the reference kernel.

    python3 tools/feyncalc/bench_fc.py                    # both kernels
    python3 tools/feyncalc/bench_fc.py --ref-cache F.tsv  # reuse a saved reference run
    python3 tools/feyncalc/bench_fc.py --binary ./other   # A/B a second Mathilda build

Runs tools/feyncalc/bench_fc.m in each kernel and prints, per case, the cold
(first call) and warm (best repeat) times and the speed-up T_ref / T_mine
(> 1: Mathilda is faster), plus the geometric mean of the warm and of the
cold speed-ups.  A LeafCount mismatch is flagged: a fast wrong answer is not a
speed-up.  The reference run is saved (--save-ref) so an A/B loop on Mathilda
alone takes seconds."""
import argparse
import math
import os
import shutil
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BENCH = os.path.join(REPO, "tools", "feyncalc", "bench_fc.m")
# The reference kernel: --ref, else $WOLFRAMSCRIPT, else wolframscript on PATH.
WOLFRAM = os.environ.get("WOLFRAMSCRIPT") or shutil.which("wolframscript") or ""


def run(binary, reps):
    env = dict(os.environ, FC_ROOT=os.path.join(REPO, "packages", "feyncalc"),
               FC_BENCH_REPS=str(reps))
    out = subprocess.run([binary, "-file", BENCH], cwd=REPO, env=env,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.decode("utf-8", "replace")
    return parse(out)


def parse(text):
    res = {}
    for line in text.splitlines():
        f = line.split("\t")
        if f[0] == "BENCH" and len(f) >= 4:
            res[f[1]] = (int(f[2]) / 1e6, int(f[3]) / 1e6, f[4] if len(f) > 4 else "")
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--binary", default=os.path.join(REPO, "Mathilda"))
    ap.add_argument("--ref", default=WOLFRAM, help="reference kernel binary")
    ap.add_argument("--ref-cache", help="read the reference run from this file instead of running it")
    ap.add_argument("--save-ref", help="save the raw reference output here")
    ap.add_argument("--reps", type=int, default=3)
    a = ap.parse_args()

    if a.ref_cache and os.path.isfile(a.ref_cache):
        ref = parse(open(a.ref_cache, encoding="utf-8").read())
    else:
        if not a.ref:
            print("bench_fc: no reference kernel (pass --ref, set WOLFRAMSCRIPT, or --ref-cache)")
            return 2
        ref = run(a.ref, a.reps)
        if a.save_ref:
            with open(a.save_ref, "w", encoding="utf-8") as fh:
                for k, (c, w, lc) in ref.items():
                    fh.write(f"BENCH\t{k}\t{round(c * 1e6)}\t{round(w * 1e6)}\t{lc}\n")
    mine = run(a.binary, a.reps)

    print(f"{'case':<16}{'mine cold':>11}{'ref cold':>11}{'x':>8}{'mine warm':>11}{'ref warm':>11}{'x':>8}  check")
    lw, lc = [], []
    for k in mine:
        if k not in ref:
            continue
        mc, mw, ml = mine[k]
        rc, rw, rl = ref[k]
        xc = max(rc, 1e-5) / max(mc, 1e-5)
        xw = max(rw, 1e-5) / max(mw, 1e-5) if k != "load" else float("nan")
        flag = "" if ml == rl else f"LeafCount {ml} vs {rl}"
        if k != "load":
            lc.append(math.log(xc))
            lw.append(math.log(xw))
        print(f"{k:<16}{mc:>11.4f}{rc:>11.4f}{xc:>8.2f}{mw:>11.4f}{rw:>11.4f}{xw:>8.2f}  {flag}")
    if lw:
        print(f"\ngeomean speed-up (excluding load): cold {math.exp(sum(lc) / len(lc)):.3f}x, "
              f"warm {math.exp(sum(lw) / len(lw)):.3f}x over {len(lw)} cases")
    missing = sorted(set(ref) - set(mine))
    if missing:
        print("missing in this kernel (crashed or erred):", ", ".join(missing))
    return 0


if __name__ == "__main__":
    sys.exit(main())
