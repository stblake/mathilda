#!/usr/bin/env python3
"""Experiment 4 -- Auto-compilation: the compiler runs without Compile[].

Runs the same workloads as ``auto_compilation.m``.  See ``README.md``.

    python3 auto_compilation.py

WHAT THIS COLUMN IS.  The experiment's subject -- a CAS noticing that ordinary
user code is numeric and compiling it without being asked -- has no Python
counterpart, because SciPy's quadrature and ODE routines are already compiled
C and never saw a Python-level expression to optimise.  So these rows are not
"the same feature, three ways"; they are the standard library each system
would actually be measured against.

That makes two rows worth reading closely.  ``solve_ivp`` is a Python callback
per step, which is why the CAS beat it by a wide margin once they compile the
right-hand side.  ``CubicSpline`` is entirely compiled C and is correspondingly
hard to beat.
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

import math

from scipy import integrate, interpolate, optimize


def f(x):
    return math.sin(x) * math.exp(-x / 10.0) + math.sqrt(x + 1.0)


def fv(x):
    return np.sin(x) * np.exp(-x / 10.0) + np.sqrt(x + 1.0)


def lorenz(t, s):
    return [10.0 * (s[1] - s[0]),
            s[0] * (28.0 - s[2]) - s[1],
            s[0] * s[1] - (8.0 / 3.0) * s[2]]


def main():
    print("Experiment 4 -- auto-compilation")
    print("")
    print("-- numeric routines over a user function --")

    bench("quad(f, 0, 50)", lambda: integrate.quad(f, 0, 50))
    check("quad(f, 0, 50)", round(integrate.quad(f, 0, 50)[0], 6))

    ks = np.arange(1, 2001, dtype=float)
    bench("sum(f(k)/k^2, k = 1..2000)", lambda: (fv(ks) / ks ** 2).sum())
    bench("prod(1 + 1/k^2, k = 1..2000)", lambda: np.prod(1.0 + 1.0 / ks ** 2))
    bench("brentq(f(x) - 2)", lambda: optimize.brentq(lambda x: f(x) - 2, 0.1, 5))

    print("")
    print("-- elementwise and scan over a machine-number list --")
    lst = np.random.rand(10 ** 6)

    def nest_scalar(m):
        x = 0.31
        for _ in range(m):
            x = 3.5 * x * (1.0 - x)
        return x

    bench("nest 3.5 x (1-x), 10^6 (CPython scalar loop)",
          lambda: nest_scalar(10 ** 6), reps=1)
    bench("cumsum(sin(list))", lambda: np.cumsum(np.sin(lst)))
    bench("sin(list)**2 + 1", lambda: np.sin(lst) ** 2 + 1.0)
    bench("cumsum(list)", lambda: np.cumsum(lst))

    print("")
    print("-- an ODE and an interpolation --")
    bench("solve_ivp, Lorenz to t = 200",
          lambda: integrate.solve_ivp(lorenz, (0.0, 200.0), [1.0, 1.0, 1.0],
                                      rtol=1e-8, atol=1e-8), reps=1)

    ni = 10 ** 4
    xs = np.linspace(0.0, 100.0, ni)
    ifn = interpolate.CubicSpline(xs, np.sin(xs))
    pts = np.random.rand(ni) * 100.0
    bench("CubicSpline, 10^4 nodes, 10^4 evaluations", lambda: ifn(pts))
    check("CubicSpline at 3.3", round(float(ifn(3.3)), 6))


if __name__ == "__main__":
    main()
