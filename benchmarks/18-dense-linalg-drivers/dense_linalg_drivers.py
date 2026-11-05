#!/usr/bin/env python3
"""Experiment 18 -- Dense linear algebra drivers (scipy.linalg column).

Same seven kernels as ``dense_linalg_drivers.m``, same order and sizes.

On this host numpy links the SAME Accelerate BLAS Mathilda does, so on these rows
all three systems call byte-identical kernels and any spread is PURE OVERHEAD --
the cleanest measurement in the suite.

Checks use the exactly reproducible SPD integer matrix from data.py, so the three
systems compare the same numbers rather than three different random matrices.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np

from harness import bench, check, require, seed
from data import spd_int_matrix

require(["numpy.linalg:solve", "numpy.linalg:inv", "numpy.linalg:det",
         "numpy.linalg:eigvals", "numpy.linalg:svd", "numpy.linalg:qr",
         "scipy.linalg:lu"])

seed()
n = 600
a = np.random.random((n, n))
b = np.random.random((n, n))
v = np.random.random(n)
s8 = spd_int_matrix(8)


def r(x, p):
    return int(math.floor(float(x) * 10 ** p + 0.5))


bench("Dot 600x600", lambda: a @ b)
check("Dot 600x600", r((s8 @ s8).sum(), 6))

bench("LinearSolve 600x600", lambda: np.linalg.solve(a, v))
check("LinearSolve 600x600",
      r(np.linalg.solve(s8, np.arange(1.0, 9.0)).sum(), 6))

bench("Inverse 600x600", lambda: np.linalg.inv(a))
check("Inverse 600x600", r(np.linalg.inv(s8).sum(), 6))

bench("Det 600x600", lambda: np.linalg.det(a))
check("Det 600x600", r(np.linalg.det(s8), 4))

bench("Eigenvalues 600x600 (general)", lambda: np.linalg.eigvals(a), reps=1)
check("Eigenvalues 600x600 (general)",
      r(np.sort(np.linalg.eigvals(s8).real).sum(), 4))

m400 = np.random.random((400, 400))
bench("SingularValueDecomposition 400x400", lambda: np.linalg.svd(m400), reps=1)
# Wolfram's SVD part 2 is the diagonal sigma MATRIX; its total is the sum of the
# singular values.
check("SingularValueDecomposition 400x400",
      r(np.linalg.svd(s8, compute_uv=False).sum(), 4))

bench("QRDecomposition 400x400", lambda: np.linalg.qr(m400), reps=1)
# Q and R signs are not canonical across libraries, so the check is |sum(R)|.
check("QRDecomposition 400x400", r(abs(np.linalg.qr(s8)[1].sum()), 4))

# ---- the decompositions and matrix functions added in v0.6180-v0.6188 -------
#
# Same cases, labels and checks as the .m half.  Minors, Adjugate and
# Orthogonalize have no scipy one-liner, so the numpy equivalent is spelled out:
# Minors as one batched det over every k-subset (itertools.combinations),
# Adjugate as det * inv, Orthogonalize as Gram-Schmidt == QR of the transpose
# with R's diagonal made positive.

import itertools

import scipy.linalg as sla

c4 = np.array([[4, 1 + 1j, 0, 0.5j], [1 - 1j, 5, 2, 0], [0, 2, 6, 1 - 0.5j],
               [-0.5j, 0, 1 + 0.5j, 7]])
g4 = np.array([[3, 1, 0, 2], [2, 4, 1, 0], [0, 1, 5, 1], [1, 0, 1, 6]], dtype=float)

spd6 = a @ a.T + n * np.eye(n)
_c = np.random.random((400, 400)) + 1j * np.random.random((400, 400))
h4 = _c @ _c.conj().T + 400 * np.eye(400)
g3 = np.random.random((300, 300))
c3 = np.random.random((300, 300)) + 1j * np.random.random((300, 300))
g2 = np.random.random((200, 200))
p1 = np.random.random((100, 100)) + 100 * np.eye(100)
r1 = np.random.random((100, 100)) / 10
m10 = np.random.random((10, 10))

require(["scipy.linalg:cholesky", "scipy.linalg:hessenberg", "scipy.linalg:eig",
         "scipy.linalg:svdvals", "scipy.linalg:expm", "scipy.linalg:logm",
         "scipy.linalg:sqrtm", "scipy.linalg:cosm",
         "scipy.linalg:fractional_matrix_power"])


def minors(m, k):
    """Every k x k minor, rows and columns in lexicographic subset order."""
    idx = np.array(list(itertools.combinations(range(m.shape[0]), k)))
    return np.linalg.det(m[idx[:, None, :, None], idx[None, :, None, :]])


def adjugate(m):
    return np.linalg.det(m) * np.linalg.inv(m)


def orthogonalize(m):
    """Gram-Schmidt on the rows: Q of QR(m^T) with diag(R) > 0, transposed."""
    q, rr = np.linalg.qr(m.T)
    return (q * np.sign(np.diag(rr))).T


bench("CholeskyDecomposition 600x600", lambda: sla.cholesky(spd6))
check("CholeskyDecomposition 600x600", r(sla.cholesky(s8).sum(), 6))

bench("CholeskyDecomposition 400x400 complex Hermitian", lambda: sla.cholesky(h4))
check("CholeskyDecomposition 400x400 complex Hermitian",
      r(sla.cholesky(c4).sum().real, 6))

m400b = np.random.random((400, 400))
bench("HessenbergDecomposition 400x400",
      lambda: sla.hessenberg(m400b, calc_q=True), reps=1)
check("HessenbergDecomposition 400x400",
      r(np.abs(sla.hessenberg(s8)).sum(), 6))

bench("Eigensystem 300x300 (general)", lambda: sla.eig(g3), reps=1)
_w, _v = sla.eig(g4)
check("Eigensystem 300x300 (general)",
      "{%d, %d}" % (r(np.abs(_w).max(), 6), r((np.abs(_v) ** 2).sum(), 6)))

bench("Eigensystem 300x300 complex", lambda: sla.eig(c3), reps=1)
check("Eigensystem 300x300 complex", r(np.abs(np.linalg.eigvals(c4)).sum(), 6))

m400c = np.random.random((400, 400))
bench("SingularValueList 400x400", lambda: sla.svdvals(m400c))
check("SingularValueList 400x400", r(sla.svdvals(s8).sum(), 6))

bench("SingularValueList 300x300 complex", lambda: sla.svdvals(c3))
check("SingularValueList 300x300 complex", r(sla.svdvals(c4).sum(), 6))

bench("Minors[m, 5] 10x10", lambda: minors(m10, 5))
check("Minors[m, 5] 10x10", r(minors(s8, 4).sum(), 3))

bench("Adjugate 200x200", lambda: adjugate(g2))
check("Adjugate 200x200", r(adjugate(s8).sum(), 3))

bench("Orthogonalize 200x200", lambda: orthogonalize(g2))
check("Orthogonalize 200x200", r(orthogonalize(s8).sum(), 6))

bench("MatrixExp 300x300", lambda: sla.expm(g3 / 10))
check("MatrixExp 300x300", r(sla.expm(s8 / 10).sum(), 6))

bench("MatrixExp 300x300 complex", lambda: sla.expm(c3 / 10))
check("MatrixExp 300x300 complex", r(sla.expm(c4 / 10).sum().real, 6))

bench("MatrixLog 100x100", lambda: sla.logm(p1), reps=1)
check("MatrixLog 100x100", r(sla.logm(s8).sum().real, 6))

bench("MatrixFunction[Sqrt] 100x100", lambda: sla.sqrtm(p1), reps=1)
check("MatrixFunction[Sqrt] 100x100", r(sla.sqrtm(s8).sum().real, 6))

bench("MatrixFunction[Cos] 100x100", lambda: sla.cosm(r1), reps=1)
check("MatrixFunction[Cos] 100x100", r(sla.cosm(s8 / 10).sum(), 6))

bench("MatrixPower[m, 0.3] 100x100",
      lambda: sla.fractional_matrix_power(p1, 0.3), reps=1)
check("MatrixPower[m, 0.3] 100x100",
      "{%d, %d}" % (r(sla.fractional_matrix_power(s8, 0.3).sum().real, 6),
                    r(sla.fractional_matrix_power(g4, 0.3).sum().real, 6)))
