#!/usr/bin/env python3
"""Experiment 7 -- Reaching the vendor BLAS and LAPACK kernels

Runs the same calls as ``blas_lapack_routing.m``, in the same order, at the
same sizes.  See ``README.md``.

    python3 blas_lapack_routing.py

On this host NumPy links the same Apple Accelerate BLAS that Mathilda does, so
on the rows that reach it these three columns are running byte-identical
kernels.  Any spread is overhead, not arithmetic.
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
# dgemm -- the row every system reaches, and therefore the calibration for the
# rest.
n=1000
A=np.random.rand(n,n)
Bm=np.random.rand(n,n)

# ---- solve -- LinearSolve, 1000x1000
# dgesv.
n=1000
A=np.random.rand(n,n)
bv=np.random.rand(n)

# ---- inverse -- Inverse, 500x500
# dgetrf + dgetri.
n=500
A=np.random.rand(n,n)

# ---- det -- Det, 500x500
# An LU factorisation and a product of the diagonal.

# ---- qr -- QRDecomposition, 500x500
# dgeqrf + dorgqr.

# ---- eigen -- Eigenvalues, 300x300 symmetric
# dsyevd for a symmetric matrix. The ordering convention is the reason this
# one is hard to route: LAPACK's order is not Wolfram's, and parity must be
# preserved.
n=300
A=np.random.rand(n,n)
S=(A+A.T)/2.0

# ---- svd -- SingularValueDecomposition, 300x300
# dgesdd.
n=300
A=np.random.rand(n,n)


def main():
    print("Experiment 7 -- Reaching the vendor BLAS and LAPACK kernels")
    print("")
    bench("Matrix multiply, 1000x1000", lambda: A @ Bm)
    bench("LinearSolve, 1000x1000", lambda: np.linalg.solve(A,bv))
    bench("Inverse, 500x500", lambda: np.linalg.inv(A))
    bench("Det, 500x500", lambda: np.linalg.det(A))
    bench("QRDecomposition, 500x500", lambda: np.linalg.qr(A))
    bench("Eigenvalues, 300x300 symmetric", lambda: np.linalg.eigvalsh(S))
    bench("SingularValueDecomposition, 300x300", lambda: np.linalg.svd(A))


if __name__ == "__main__":
    main()
