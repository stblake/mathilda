#!/usr/bin/env python3
"""Experiment 3 -- Optimising codegen: CSE, folding, DCE, LICM.

Runs the same bodies as ``compile_optimising_codegen.m``.  See ``README.md``.

    python3 compile_optimising_codegen.py

WHAT THIS COLUMN IS FOR.  The experiment is about a bytecode compiler's
optimiser, which Python does not have an equivalent of -- CPython's peephole
optimiser does not do CSE or LICM across a loop.  These rows are therefore the
scalar-loop baseline the CAS optimisers are trying to beat, not a competing
optimiser.  **numba is not installed on this host**; with it, the first four
rows would be a genuine optimiser comparison.

The last row IS vectorisable and is a fair NumPy result.
"""

import time

import numpy as np


def bench(label, fn, reps=3):
    """One untimed warm-up, then the MINIMUM of `reps` timed runs."""
    fn()
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    print("%-52s%s ms" % (label, round(1000.0 * min(ts), 3)))


def check(label, value):
    print("%-52scheck = %s" % (label, value))


def horner(x):
    """A degree-40 Horner chain with constant coefficients."""
    s = 1.0
    for k in range(1, 41):
        s = s * x + (1.0 + 1.0 / (k + 1.0))
    return s


def cse(x, y):
    """A repeated subexpression: x*y appears twice."""
    import math
    return x * y + math.sin(x * y)


def licm(x, n):
    """A loop-invariant subexpression inside a while loop."""
    import math
    s, k = 0.0, 0
    while k < n:
        s = s + x * math.sqrt(2.0) + math.sin(1.0)
        k += 1
    return s


def newt(a):
    """Newton's method -- the micro-benchmark the codegen work was tuned on."""
    x, k = 1.0, 0
    while k < 20:
        x = x - (x * x - a) / (2.0 * x)
        k += 1
    return x


def main():
    print("Experiment 3 -- optimising codegen")
    print("")
    bench("Horner, degree 40, 10^5 calls (CPython)",
          lambda: [horner(0.5) for _ in range(10 ** 5)], reps=1)
    check("Horner, degree 40", horner(0.5))

    bench("x*y + sin(x*y), 10^5 calls (CPython)",
          lambda: [cse(0.5, 0.25) for _ in range(10 ** 5)], reps=1)
    check("x*y + sin(x*y)", cse(0.5, 0.25))

    bench("loop-invariant body, 10^6 iterations (CPython)",
          lambda: licm(0.5, 10 ** 6), reps=1)
    check("loop-invariant body", licm(0.5, 10))

    bench("Newton, 20 iterations, 10^5 calls (CPython)",
          lambda: [newt(2.0) for _ in range(10 ** 5)], reps=1)
    check("Newton (sqrt 2)", newt(2.0))

    x = np.linspace(0.0, 1.0, 10 ** 6)
    bench("degree-5 polynomial over 10^6 points (NumPy)",
          lambda: 1.0 + 2.0 * x + 3.0 * x ** 2 + 4.0 * x ** 3
                  + 5.0 * x ** 4 + 6.0 * x ** 5)


if __name__ == "__main__":
    main()
