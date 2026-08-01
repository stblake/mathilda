#!/usr/bin/env python3
"""Experiment 9 -- First sweep: classical HPC primitives

A representative kernel from each of the seven groups the first sweep
covered.  See ``README.md``; the full 43-kernel sweep is in
``comparisons/hpc_bench.py``, from which this file is generated.

    python3 hpc_sweep_primitives.py

Where a row has no NumPy equivalent -- arbitrary precision, symbolic rule
dispatch -- it says so rather than substituting something that is not the same
computation.
"""

import numpy as np

import time


def bench(label, fn, reps=3):
    """One untimed warm-up, then the MINIMUM of `reps` timed runs.

    The minimum, not the mean: we are measuring the cost of the work, and every
    source of noise on a loaded machine can only add.
    """
    fn()
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    print("%-52s%s ms" % (label, round(1000.0 * min(ts), 3)))


def bench1(label, fn):
    """A single timed run, no warm-up, for kernels that take seconds."""
    t0 = time.perf_counter()
    fn()
    print("%-52s%s ms  (1 run)" % (label, round(1000.0 * (time.perf_counter() - t0), 3)))


def check(label, value):
    print("%-52scheck = %s" % (label, value))


def nest(f, x, m):
    """Nest[f, x, m] -- apply f to x, m times."""
    for _ in range(m):
        x = f(x)
    return x

# ---- matmul -- Matrix multiply, 1000x1000
# Dense linear algebra: dgemm, which all three systems reach.
n=1000
A=np.random.rand(n,n)
Bm=np.random.rand(n,n)

# ---- fft -- Fourier, 2^20 reals
# Spectral: one large transform, dominated by FFTW.
v=np.random.rand(2**20)

# ---- triad -- STREAM triad, a = b + 3 c
# Array and memory: the STREAM triad, this machine's bandwidth reference.
n=10**7
x=np.random.rand(n)
y=np.random.rand(n)

# ---- sort -- Sort
# Order statistics at 10^7.

# ---- jacobi -- Jacobi 5-point relaxation, 512^2, 100 sweeps
# Stencils: a 5-point relaxation, the shape every PDE solver is built from.
n=512
u0=np.random.rand(n,n)
def jac(u):
    return (np.roll(u,-1,0)+np.roll(u,1,0)+np.roll(u,-1,1)+np.roll(u,1,1))/4.0

# ---- logistic -- Logistic map, 10^7 iterations
# A scalar kernel via Compile[]: a tight dependent loop with no array in it at
# all.
def lg(x0, m):
    x = x0
    for _ in range(m): x = 3.9*x*(1.0-x)
    return x

# ---- mandel -- Mandelbrot, 800x800, 100 iterations
# The same, over a 2-D domain with an early exit.
n=800
h=2.5/(n-1)
def mandelgrid(m):
    cy = (-1.25 + h*np.arange(n))[:,None]
    cx = (-2.0 + h*np.arange(n))[None,:]
    zx = np.zeros((n,n)); zy = np.zeros((n,n))
    it = np.zeros((n,n), dtype=np.int64)
    live = np.ones((n,n), dtype=bool)
    for _ in range(m):
        live &= (zx*zx + zy*zy) < 4.0
        if not live.any(): break
        t = zx*zx - zy*zy + cx
        zy = np.where(live, 2.0*zx*zy + cy, zy)
        zx = np.where(live, t, zx)
        it += live
    return it

# ---- sieve -- Sieve of Eratosthenes to 10^7
# Integer and combinatorial.
def sv(m):
    s = np.ones(m+1, dtype=bool); s[:2] = False
    i = 2
    while i*i <= m:
        if s[i]: s[i*i::i] = False
        i += 1
    return int(s.sum())

# ---- fib -- Naive recursive Fibonacci, fib(25)
# Rule dispatch rather than arithmetic -- the evaluator's own speed, and the
# one row where CPython beats both CAS.
def fib(k):
    return k if k < 2 else fib(k-1)+fib(k-2)

# ---- pi -- pi to 100,000 digits
# Arbitrary precision: MPFR, where both CAS call the same library and NumPy
# has no equivalent.
import mpmath as mp
def mppi(d):
    mp.mp.dps = d + 10
    return +mp.pi

# ---- fact -- 50000! (exact)
import math as _m
import sys as _s
_s.set_int_max_str_digits(0)

# ---- bigmul -- Product of two 10^6-bit integers
import sys as _s
_s.set_int_max_str_digits(0)
p=2**1000003 - 1
q=3**631305


def main():
    print("Experiment 9 -- First sweep: classical HPC primitives")
    print("")
    bench("Matrix multiply, 1000x1000", lambda: A @ Bm)
    bench("Fourier, 2^20 reals", lambda: np.fft.fft(v))
    bench("STREAM triad, a = b + 3 c", lambda: x + 3.0*y)
    bench("Sort", lambda: np.sort(x))
    bench("Jacobi 5-point relaxation, 512^2, 100 sweeps", lambda: nest(jac,u0,100))
    bench1("Logistic map, 10^7 iterations", lambda: lg(0.5, 10**7))
    check("Logistic map, 10^7 iterations", lg(0.5, 10**7))
    bench1("Mandelbrot, 800x800, 100 iterations", lambda: mandelgrid(100))
    check("Mandelbrot, 800x800, 100 iterations", int(mandelgrid(100).sum()))
    bench1("Sieve of Eratosthenes to 10^7", lambda: sv(10**7))
    check("Sieve of Eratosthenes to 10^7", sv(10**7))
    bench1("Naive recursive Fibonacci, fib(25)", lambda: fib(25))
    check("Naive recursive Fibonacci, fib(25)", fib(25))
    bench1("pi to 100,000 digits", lambda: mppi(100000 + kk))
    check("pi to 100,000 digits", int(mp.nint((mppi(100000) - 3) * mp.mpf(10)**20)))
    bench1("50000! (exact)", lambda: _m.factorial(50000 + kk))
    check("50000! (exact)", len(str(_m.factorial(50000))))
    bench1("Product of two 10^6-bit integers", lambda: ((p + kk)*q).bit_length())
    check("Product of two 10^6-bit integers", len(str(p*q)))


if __name__ == "__main__":
    main()
