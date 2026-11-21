#!/usr/bin/env python3
"""Run tests/data/tc_unwind_stress.m in parallel copies of ./Mathilda.

Every copy must exit 0 and print "3000 True".  A crash (the TimeConstrained
async-unwind class fixed in v0.6662/v0.6664) shows as a signal exit; the run is
parallel on purpose, since the failure needs timeouts to land at random points
and load spreads them.  Exit 77 (skip) when ./Mathilda is not built.
"""
import os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MATH = os.path.join(ROOT, "Mathilda")
SCRIPT = os.path.join(ROOT, "tests", "data", "tc_unwind_stress.m")
COPIES = int(os.environ.get("TC_UNWIND_COPIES", "12"))

def one(i):
    try:
        p = subprocess.run([MATH, "-file", SCRIPT], cwd=ROOT, capture_output=True, text=True, timeout=900)
        return i, p.returncode, p.stdout.strip().splitlines()[-1:] or [""]
    except subprocess.TimeoutExpired:
        return i, "timeout", [""]

def main():
    if not os.access(MATH, os.X_OK):
        print("SKIP: ./Mathilda not built"); sys.exit(77)
    with ThreadPoolExecutor(max_workers=COPIES) as ex:
        res = list(ex.map(one, range(COPIES)))
    bad = [(i, rc, out[0]) for i, rc, out in res if rc != 0 or out[0] != "3000 True"]
    for i, rc, out in bad:
        print(f"copy {i}: exit {rc}, last line {out!r}")
    print(f"{COPIES - len(bad)}/{COPIES} copies clean")
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
