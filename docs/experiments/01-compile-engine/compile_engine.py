#!/usr/bin/env python3
"""Experiment 1 -- Compile[]: a typed bytecode VM for numeric code.

Runs the same four kernels as ``compile_engine.m``.  See ``README.md``.

    python3 compile_engine.py

WHERE THIS COLUMN IS NOT A LIBRARY RESULT, AND WHY IT SAYS SO.  Three of the
four kernels are tight SCALAR loops.  Python's answer to those is numba or
Cython, and **numba is not installed on this host**, so the first three rows
are plain CPython and must be read as "an interpreted scalar loop", not as
"what Python can do".  They are the honest comparison for what the two CAS are
doing on the same rows -- both of which are also starting from an interpreter
and compiling.

The fourth kernel IS vectorisable, and its NumPy column is a fair library
result.
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


def nest(f, x, m):
    for _ in range(m):
        x = f(x)
    return x


def logi(n):
    """A scalar recurrence: nothing to vectorise, each iterate needs the last."""
    x, k = 0.31, 0
    while k < n:
        x = 3.9 * x * (1.0 - x)
        k += 1
    return x


def mand(w, maxit):
    """A 2-D domain with a data-dependent early exit."""
    c = 0
    for i in range(w):
        cr = -2.0 + 3.0 * i / w
        for j in range(w):
            ci = -1.5 + 3.0 * j / w
            zr = zi = 0.0
            k = 0
            while k < maxit and zr * zr + zi * zi < 4.0:
                zr, zi = zr * zr - zi * zi + cr, 2.0 * zr * zi + ci
                k += 1
            c += k
    return c


def lj(n):
    """An all-pairs scalar sum."""
    import math
    e = 0.0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            dx = math.sin(1.0 * i) - math.sin(1.0 * j)
            dy = math.cos(1.3 * i) - math.cos(1.3 * j)
            dz = math.sin(0.7 * i) - math.sin(0.7 * j)
            r2 = dx * dx + dy * dy + dz * dz + 0.5
            ir6 = 1.0 / (r2 * r2 * r2)
            e += 4.0 * (ir6 * ir6 - ir6)
    return e


def mcpi(m):
    """The one vectorisable kernel of the four -- a fair NumPy result."""
    u = np.random.rand(m)
    v = np.random.rand(m)
    return 4.0 * float((u * u + v * v <= 1.0).sum()) / m


def main():
    print("Experiment 1 -- Compile[]: a typed bytecode VM")
    print("")
    bench("logistic map, 10^7 iterations (CPython)", lambda: logi(10 ** 7), reps=1)
    check("logistic map, 10^7 iterations", logi(10 ** 7))
    bench("Mandelbrot, 800^2, 100 iterations (CPython)",
          lambda: mand(800, 100), reps=1)
    check("Mandelbrot, 800^2, 100 iterations", mand(800, 100))
    bench("Lennard-Jones energy, 1452 bodies (CPython)", lambda: lj(1452), reps=1)
    check("Lennard-Jones energy, 1452 bodies", lj(1452))
    bench("Monte Carlo pi, 10^7 samples (vectorised NumPy)",
          lambda: mcpi(10 ** 7))


if __name__ == "__main__":
    main()
