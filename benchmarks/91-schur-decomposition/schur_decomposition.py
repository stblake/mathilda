#!/usr/bin/env python3
"""Experiment 30 -- Schur decomposition (scipy.linalg column).

Same cases, order, and sizes as schur_decomposition.m.  scipy.linalg.schur wraps
the same Accelerate dgees/zgees Mathilda calls, and scipy.linalg.qz the same
dgges/zgges, so on this host the kernels are byte-identical and any spread is
pure overhead.  The "packed" row is scipy's identical schur() -- it exists to
compare against Mathilda's packed/NDArray fast path, which no longer round-trips
the result through boxed Exprs.

Checks reconstruct a reproducible symmetric matrix; residual is 0 for any correct
(non-unique) Schur basis, so the two systems agree without pinning a basis.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np

from scipy.linalg import schur, qz

from harness import bench, check, require, seed
from data import spd_int_matrix

require(["scipy.linalg:schur", "scipy.linalg:qz"])

seed()


def r(x, p):
    return int(math.floor(float(x) * 10 ** p + 0.5))


s8 = np.asarray(spd_int_matrix(8), dtype=float)
s8b = s8 + np.eye(8)
c8 = s8 + 1j * s8b


def schur_resid(m):
    out = "complex" if np.iscomplexobj(m) else "real"
    T, Z = schur(m, output=out)
    return r(np.max(np.abs(m - Z @ T @ Z.conj().T)), 6)


def gen_resid(m, a):
    S, T, Q, Zg = qz(m, a)
    return r(max(np.max(np.abs(m - Q @ S @ Zg.conj().T)),
                 np.max(np.abs(a - Q @ T @ Zg.conj().T))), 6)


n = 300
bench("SchurDecomposition 300x300", lambda: schur(np.random.random((n, n))), reps=1)
check("SchurDecomposition 300x300", schur_resid(s8))

n2 = 200
bench("SchurDecomposition complex 200x200",
      lambda: schur(np.random.random((n2, n2)) + 1j * np.random.random((n2, n2)),
                    output="complex"), reps=1)
check("SchurDecomposition complex 200x200", schur_resid(c8))

bench("SchurDecomposition generalized 200x200",
      lambda: qz(np.random.random((n2, n2)), np.random.random((n2, n2))), reps=1)
check("SchurDecomposition generalized 200x200", gen_resid(s8, s8b))
