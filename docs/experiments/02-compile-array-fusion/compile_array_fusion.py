#!/usr/bin/env python3
"""Experiment 2 -- Compiled machine arrays and fused elementwise loops.

Runs the same expressions as ``compile_array_fusion.m``.  See ``README.md``.

    python3 compile_array_fusion.py

WHAT THE NUMPY COLUMN IS HERE.  NumPy is the UNFUSED reference: ``a**2 + 2*a +
1`` in NumPy is four passes over 8 MB with three temporaries, exactly like the
interpreted CAS form, because NumPy has no fusion either.  That is what makes
it the right control for this experiment -- the compiled-and-fused CAS numbers
should beat it, and the point of the row is by how much.

(``numexpr`` and ``numba`` do fuse; neither is installed on this host, and
substituting one would be comparing a different library.)
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

N = 10 ** 6
v = np.random.rand(N)


def main():
    print("Experiment 2 -- compiled machine arrays and fused elementwise loops")
    print("")
    print("-- NumPy is UNFUSED: each of these is several passes --")
    bench("v**2 + 2*v + 1", lambda: v ** 2 + 2 * v + 1)
    bench("(v**2 + 2*v + 1).sum()", lambda: (v ** 2 + 2 * v + 1).sum())
    bench("(v + v*v).sum()", lambda: (v + v * v).sum())
    bench("sqrt(v) + v**2", lambda: np.sqrt(v) + v ** 2)

    print("")
    print("-- libm-bound bodies --")
    bench("sin(v)*exp(-v) + sqrt(v)",
          lambda: np.sin(v) * np.exp(-v) + np.sqrt(v))
    bench("(sin(v)*exp(-v) + sqrt(v)).sum()",
          lambda: (np.sin(v) * np.exp(-v) + np.sqrt(v)).sum())

    print("")
    print("-- the elementwise primitives, for scale --")
    bench("sin(v)", lambda: np.sin(v))
    bench("exp(v)", lambda: np.exp(v))
    bench("v + 3.0*v", lambda: v + 3.0 * v)

    check("(v**2 + 2*v + 1).sum()", round(float((v ** 2 + 2 * v + 1).sum()), 6))


if __name__ == "__main__":
    main()
