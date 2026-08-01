#!/usr/bin/env python3
"""Experiment 5 -- Machine integers: int64 as a peer of float64

Runs the same algorithms as ``machine_integers.m``, in the same order, with
the same sizes.  See ``README.md``.

    python3 machine_integers.py

WHERE THE PYTHON COLUMN IS NOT NUMPY, AND WHY IT SAYS SO.  The sieve and the
Game of Life vectorise, so those are a fair library comparison.  Collatz and
the recursive Fibonacci do NOT vectorise -- one is a data-dependent scalar
loop, the other is recursion -- so their Python column is plain CPython.
**numba is not installed on this host**, so those two rows must be read as
"an interpreted scalar loop", not as "what Python can do".  The recursive
Fibonacci is included as the control that separates array performance from
evaluator performance.
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

# ---- sieve -- Sieve of Eratosthenes to 10^7
# A pure integer-array workload inside Compile[]: a flag vector, strided
# writes, and a count. The compiled form is the one that was measured; an
# interpreted strided Part assignment is a different and far slower program.
def sv(m):
    s = np.ones(m+1, dtype=bool); s[:2] = False
    i = 2
    while i*i <= m:
        if s[i]: s[i*i::i] = False
        i += 1
    return int(s.sum())

# ---- collatz -- Collatz longest chain below 10^6
# A data-dependent scalar loop -- the case Compile[] exists for, and the one
# where an integer type must not silently become a float: 3n+1 leaves a
# double's exact range long before it leaves int64.
def cz(m):
    bl = 0
    for k in range(1, m+1):
        q = k; ln = 0
        while q > 1:
            q = q//2 if q % 2 == 0 else 3*q+1
            ln += 1
        if ln > bl: bl = ln
    return bl

# ---- life -- Game of Life, 256^2, 100 generations
# Eight rotations, an integer sum and two comparisons. The neighbour count is
# an int64 grid, so its UnitStep sees integers and not reals -- which is what
# the narrowing kernel category of this experiment is for.
n=256
g0=np.random.randint(0,2,(n,n))
def nb(g):
    s = -g
    for i in (-1,0,1):
        for j in (-1,0,1):
            s = s + np.roll(np.roll(g,-i,0),-j,1)
    return s
def life(g):
    k = nb(g)
    return (k==3).astype(np.int64) + (k==2).astype(np.int64)*g

# ---- fib -- Naive recursive Fibonacci, fib(25)
# Not an array kernel at all: pure rule dispatch, and the control that says
# how much of the above is the buffer and how much is the evaluator.
def fib(k):
    return k if k < 2 else fib(k-1)+fib(k-2)


def main():
    print("Experiment 5 -- Machine integers: int64 as a peer of float64")
    print("")
    bench1("Sieve of Eratosthenes to 10^7", lambda: sv(10**7))
    check("Sieve of Eratosthenes to 10^7", sv(10**7))
    bench1("Collatz longest chain below 10^6", lambda: cz(10**6))
    check("Collatz longest chain below 10^6", cz(10**6))
    bench("Game of Life, 256^2, 100 generations", lambda: nest(life,g0,100))
    bench1("Naive recursive Fibonacci, fib(25)", lambda: fib(25))
    check("Naive recursive Fibonacci, fib(25)", fib(25))


if __name__ == "__main__":
    main()
