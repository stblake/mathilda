#!/usr/bin/env python3
"""Experiment 54 -- Eigenproblems (scipy column).

Same kernels as ``eigenproblems.m``, same order and labels.

EXECUTION comparison against compiled LAPACK (scipy.linalg) and ARPACK
(scipy.sparse.linalg).  Timing on large random matrices; checks on the small
reproducible spd_int_matrix(8) or a fixed matrix, as order-invariant scalars.

Eigensystem was absent in Mathilda until v0.6185 (the .m side benchIf-SKIPs
it on an older binary).
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math
import numpy as np
import scipy.linalg as sla
import scipy.sparse.linalg as ssla

from harness import bench, check, require, seed
from data import spd_int_matrix

require(["numpy.linalg:eigvals", "numpy.linalg:eigvalsh", "scipy.linalg:eigh",
         "scipy.linalg:eig", "scipy.sparse.linalg:eigsh"])


def r(x, p):
    return int(math.floor(float(x) * 10 ** p + 0.5))


seed()
ns = 300
m = np.random.random((ns, ns)); sy = (m + m.T) / 2               # symmetric
ge = np.random.random((ns, ns))                                  # general
mv = np.random.random((250, 250)); sv = (mv + mv.T) / 2          # eigenvectors
mb = np.random.random((500, 500)); big = (mb + mb.T) / 2         # Arnoldi

sm = spd_int_matrix(8)
gen4 = np.array([[3, 1, 0, 2], [2, 4, 1, 0], [0, 1, 5, 1], [1, 0, 1, 6]], dtype=float)
A5 = spd_int_matrix(5)
B5 = np.diag([2., 3., 4., 5., 6.])

bench("Eigenvalues symmetric 300x300", lambda: np.linalg.eigvalsh(sy), reps=1)
check("Eigenvalues symmetric 300x300", r(np.linalg.eigvalsh(sm).max(), 6))

bench("Eigenvalues general 300x300", lambda: np.linalg.eigvals(ge), reps=1)
check("Eigenvalues general 300x300", r(np.abs(np.linalg.eigvals(gen4)).sum(), 4))

bench("Eigenvectors symmetric 250x250", lambda: np.linalg.eigh(sv), reps=1)  # no check

# Generalized symmetric-definite eigenvalues -- scipy.linalg.eigh(A, B), one LAPACK
# call, vs Mathilda's symbolic char-poly path. Same 5x5 loop shape.
bench("Generalized eigenvalues 5x5 x50", lambda: [sla.eigh(A5, B5) for _ in range(50)], reps=1)
check("Generalized eigenvalues 5x5 x50", r(sla.eigh(A5, B5, eigvals_only=True).sum(), 6))

bench("Eigenvalues Arnoldi k=6 of 500x500", lambda: ssla.eigsh(big, k=6))
check("Eigenvalues Arnoldi k=6 of 500x500", r(np.linalg.eigvalsh(sm).sum(), 6))

# Eigensystem: scipy returns (w, v) together.
bench("Eigensystem symmetric 250x250", lambda: np.linalg.eigh(sv), reps=1)
check("Eigensystem symmetric 250x250", r(np.linalg.eigvalsh(sm).max(), 6))

# ---- Eigensystem (v0.6185) and the k forms ----------------------------------
# Dense LAPACK for the full Hermitian problem; ARPACK (eigsh / svds) for the
# count forms, which is what the k argument asks for.  Checks as in the .m.

c4 = np.array([[4, 1 + 1j, 0, 0.5j], [1 - 1j, 5, 2, 0], [0, 2, 6, 1 - 0.5j],
               [-0.5j, 0, 1 + 0.5j, 7]])
_c = np.random.random((250, 250)) + 1j * np.random.random((250, 250))
hc = (_c + _c.conj().T) / 2


def es_chk(w, v):
    return "{%d, %d}" % (r(np.real(w).sum(), 6), r((np.abs(v) ** 2).sum(), 6))


_sw = np.linalg.eigvalsh(sm)
_by_mag = _sw[np.argsort(-np.abs(_sw))]

bench("Eigensystem Hermitian 250x250 complex", lambda: np.linalg.eigh(hc), reps=1)
check("Eigensystem Hermitian 250x250 complex", es_chk(*np.linalg.eigh(c4)))

bench("Eigenvalues k=6 (largest) of symmetric 300x300",
      lambda: ssla.eigsh(sy, k=6, which="LM", return_eigenvectors=False), reps=1)
check("Eigenvalues k=6 (largest) of symmetric 300x300", r(_by_mag[:3].sum(), 6))

bench("Eigenvalues k=-6 (smallest) of symmetric 300x300",
      lambda: ssla.eigsh(sy, k=6, sigma=0, which="LM", return_eigenvectors=False), reps=1)
check("Eigenvalues k=-6 (smallest) of symmetric 300x300", r(_by_mag[-3:].sum(), 6))

bench("Eigensystem k=6 of symmetric 300x300", lambda: ssla.eigsh(sy, k=6, which="LM"), reps=1)
_w, _v = np.linalg.eigh(sm)
_top = np.argsort(-np.abs(_w))[:3]
check("Eigensystem k=6 of symmetric 300x300", es_chk(_w[_top], _v[:, _top]))

bench("SingularValueList k=6 of 300x300",
      lambda: ssla.svds(ge, k=6, return_singular_vectors=False), reps=1)
check("SingularValueList k=6 of 300x300", r(np.sort(sla.svdvals(sm))[-3:].sum(), 6))
