#!/usr/bin/env python3
"""Experiment 124 -- sparse linear algebra at scale (scipy column).

Same cases, labels and checks as ``sparse_linalg.m``; see its header for the
problems.  Baselines are scipy.sparse.linalg:

* direct     ``spsolve`` (SuperLU, COLAMD) -- scipy ships neither UMFPACK nor
             CHOLMOD, so SuperLU is the direct row for SPD, nonsymmetric and
             complex systems alike.
* Krylov     ``cg`` / ``gmres`` (restart 50) / ``bicgstab`` at the same
             relative tolerance.  scipy has no ILU0 or ILUT proper; its only
             incomplete factorisation is SuperLU's threshold ILU (``spilu``),
             used in natural order without pivoting: drop 1e-2 / fill 1 for
             ILU0, drop 1e-4 / fill 5 for ILUT.
* many-RHS   ``splu(A)`` once, then ``.solve`` per right-hand side.
* heat       ``expm_multiply`` (Al-Mohy-Higham), the same algorithm.
* eigen      ARPACK ``eigsh`` (``sigma=0`` shift-invert for the smallest).
* exact      sympy's ``DomainMatrix`` over ZZ/QQ (fraction-free / Bareiss,
             RREF), the exact engine behind ``Matrix``.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from harness import bench, bench_once, check, require

require(["scipy.sparse.linalg:spsolve", "scipy.sparse.linalg:splu",
         "scipy.sparse.linalg:spilu", "scipy.sparse.linalg:cg",
         "scipy.sparse.linalg:gmres", "scipy.sparse.linalg:bicgstab",
         "scipy.sparse.linalg:expm_multiply", "scipy.sparse.linalg:eigsh",
         "sympy.polys.matrices:DomainMatrix"])

res = {}


def tri(k):
    return sp.diags_array([-1.0, 2.0, -1.0], offsets=[-1, 0, 1], shape=(k, k), format="csr")


def eye(k):
    return sp.identity(k, format="csr")


def lap2(k):
    return (sp.kron(tri(k), eye(k)) + sp.kron(eye(k), tri(k))).tocsr()


def lap3(k):
    return (sp.kron(tri(k), eye(k * k)) + sp.kron(sp.kron(eye(k), tri(k)), eye(k))
            + sp.kron(eye(k * k), tri(k))).tocsr()


def rhs(n):
    return (np.arange(1, n + 1) % 7) - 3.0


def nrm(x, scale):
    return int(round(scale * float(np.linalg.norm(x))))


def spilu_op(a, kind):
    kw = dict(drop_tol=1e-2, fill_factor=1) if kind == "ILU0" else dict(drop_tol=1e-4, fill_factor=5)
    P = spla.spilu(a.tocsc(), permc_spec="NATURAL", diag_pivot_thresh=0.0, **kw)
    return spla.LinearOperator(a.shape, P.solve)


def krylov(a, b, meth, kind, tol):
    extra = {"restart": 50} if meth == "gmres" else {}
    M = spilu_op(a, kind) if kind else None
    return getattr(spla, meth)(a, b, rtol=tol, M=M, maxiter=5000, **extra)[0]


def keep(key, fn):
    return lambda: res.__setitem__(key, fn())


# ---- Poisson -----------------------------------------------------------------

L500 = lap2(500).tocsc(); b500 = rhs(500 ** 2)
bench("poisson2d 500^2: LinearSolve direct", keep("x", lambda: spla.spsolve(L500, b500)), reps=1)
check("poisson2d 500^2: LinearSolve direct", nrm(res["x"], 10))

L30 = lap3(30).tocsc(); b30 = rhs(30 ** 3)
bench_once("poisson3d 30^3: LinearSolve direct", keep("x", lambda: spla.spsolve(L30, b30)))
check("poisson3d 30^3: LinearSolve direct", nrm(res["x"], 10))

L100 = lap3(100); b100 = rhs(100 ** 3)
bench_once("poisson3d 100^3 (10^6): CG 1e-6",
           keep("x", lambda: krylov(L100, b100, "cg", None, 1e-6)))
check("poisson3d 100^3 (10^6): CG 1e-6", nrm(res["x"], 1))

# ---- convection-diffusion (nonsymmetric) ---------------------------------------

skew = sp.diags_array([-1.0, 1.0], offsets=[-1, 1], shape=(300, 300), format="csr")
C300 = (lap2(300) + 0.4 * sp.kron(eye(300), skew) + 0.2 * sp.kron(skew, eye(300))).tocsr()
b300 = rhs(300 ** 2)
bench("convdiff 300^2: LinearSolve Multifrontal",
      keep("x", lambda: spla.spsolve(C300.tocsc(), b300)), reps=1)
check("convdiff 300^2: LinearSolve Multifrontal", nrm(res["x"], 1))
bench("convdiff 300^2: GMRES + ILUT 1e-8",
      keep("x", lambda: krylov(C300, b300, "gmres", "ILUT", 1e-8)), reps=1)
check("convdiff 300^2: GMRES + ILUT 1e-8", nrm(res["x"], 1))
bench("convdiff 300^2: BiCGSTAB + ILU0 1e-8",
      keep("x", lambda: krylov(C300, b300, "bicgstab", "ILU0", 1e-8)), reps=1)
check("convdiff 300^2: BiCGSTAB + ILU0 1e-8", nrm(res["x"], 1))

# ---- complex Helmholtz-like -------------------------------------------------------

H300 = (lap2(300) + (-0.3 + 0.05j) * eye(300 ** 2)).tocsc()
bench("helmholtz 300^2 complex: LinearSolve direct",
      keep("x", lambda: spla.spsolve(H300, b300)), reps=1)
check("helmholtz 300^2 complex: LinearSolve direct", nrm(res["x"], 10))

# ---- many right-hand sides ------------------------------------------------------

L300 = lap2(300).tocsc()


def lsf100():
    lu = spla.splu(L300)
    for k in range(1, 101):
        res["x"] = lu.solve(b300 + k)
    res["lu"] = lu


bench("lsf 300^2: LinearSolveFunction, factor + 100 solves", lsf100, reps=1)
check("lsf 300^2: LinearSolveFunction, factor + 100 solves",
      nrm(res["lu"].solve(b300 + 100.0), 1))

# ---- heat step ----------------------------------------------------------------------

u500 = np.sin(np.arange(1, 500 ** 2 + 1) / 100.0)
mL500 = (-L500).tocsr()
bench("heat 500^2: MatrixExp[-L, u0] (action)", keep("y", lambda: spla.expm_multiply(mL500, u500)))
check("heat 500^2: MatrixExp[-L, u0] (action)", nrm(res["y"], 1000))

# ---- exact small matrices -------------------------------------------------------------

from sympy import QQ, ZZ, Integer, Matrix, Rational, symbols
from sympy.polys.matrices import DomainMatrix


def ex_e(n, sd, i, j):
    s = (48271 * (n * (i - 1) + j) + sd) % 2147483647
    return ((s % 19) - 9 if s % 3 == 0 else 0) + (10 if i == j else 0)


def ex_rows(n, sd, rows, cols):
    return [[ex_e(n, sd, i, j) for j in range(1, cols + 1)] for i in range(1, rows + 1)]


def dm(rows, dom=ZZ):
    return DomainMatrix([[dom(v) for v in r] for r in rows], (len(rows), len(rows[0])), dom)


def ex_r(n, r):
    return dm(ex_rows(n, 7, n, n)) * dm(ex_rows(n, 11, n, r)) * dm(ex_rows(n, 13, r, n))


m40 = dm(ex_rows(40, 1, 40, 40)); m30 = dm(ex_rows(30, 2, 30, 30))
m20 = dm(ex_rows(20, 3, 20, 20)); m10 = dm(ex_rows(10, 4, 10, 10))
r40 = ex_r(40, 30); r30 = ex_r(30, 24)
b40e = DomainMatrix([[QQ(i)] for i in range(1, 41)], (40, 1), QQ)


def qsum(m):
    """Sum of a QQ/ZZ DomainMatrix's entries, printed as p/q like InputForm."""
    t = sum((v for row in m.tolist() for v in row), Rational(0))
    return str(Rational(t))


bench("exact 40x40: Det", keep("d", lambda: m40.det()))
check("exact 40x40: Det", int(res["d"]))
bench("exact 30x30: Inverse", keep("iv", lambda: m30.convert_to(QQ).inv()))
check("exact 30x30: Inverse", qsum(res["iv"].to_Matrix()))
bench("exact 40x40: LinearSolve", keep("x", lambda: m40.convert_to(QQ).lu_solve(b40e)))
check("exact 40x40: LinearSolve", qsum(res["x"].to_Matrix()))
bench("exact 40x40: MatrixRank (rank 30)", keep("rk", lambda: r40.convert_to(QQ).rank()))
check("exact 40x40: MatrixRank (rank 30)", res["rk"])
bench("exact 30x30: NullSpace (nullity 6)", keep("ns", lambda: r30.convert_to(QQ).nullspace()))
check("exact 30x30: NullSpace (nullity 6)", res["ns"].shape[0])
bench("exact 20x20: CharacteristicPolynomial", keep("cp", lambda: m20.charpoly()))
# charpoly() is det(x I - m), highest degree first; n = 20 is even.
check("exact 20x20: CharacteristicPolynomial",
      int(sum(int(c) * 3 ** (len(res["cp"]) - 1 - k) for k, c in enumerate(res["cp"]))))
m10q = m10.convert_to(QQ)
bench("exact 10x10: Inverse x100", lambda: [res.__setitem__("iv", m10q.inv()) for _ in range(100)])
check("exact 10x10: Inverse x100", qsum(res["iv"].to_Matrix()))

# ---- eigenvalues with a count (ARPACK) ------------------------------------------------
# H = L + V, V(i, j) = 0.001 (i^2 + 2 j^2): see the .m half.

_i = np.arange(1, 201)[:, None]
_j = np.arange(1, 201)[None, :]
H200 = (lap2(200) + sp.diags_array(0.001 * (_i ** 2 + 2 * _j ** 2).ravel().astype(float))).tocsr()
bench("eigen 200^2 L + V: Eigenvalues[H, 6] (largest)",
      keep("ev", lambda: spla.eigsh(H200, k=6, which="LM", return_eigenvectors=False)), reps=1)
check("eigen 200^2 L + V: Eigenvalues[H, 6] (largest)", int(round(1e6 * float(res["ev"].sum()))))
bench("eigen 200^2 L + V: Eigenvalues[H, -6] (smallest)",
      keep("ev", lambda: spla.eigsh(H200.tocsc(), k=6, sigma=0, which="LM",
                                    return_eigenvectors=False)), reps=1)
check("eigen 200^2 L + V: Eigenvalues[H, -6] (smallest)", int(round(1e6 * float(res["ev"].sum()))))
bench("eigen 200^2 L + V: Eigensystem[H, -6]",
      keep("es", lambda: spla.eigsh(H200.tocsc(), k=6, sigma=0, which="LM")), reps=1)
check("eigen 200^2 L + V: Eigensystem[H, -6]",
      "{%d, %d}" % (round(1e6 * res["es"][0].sum()), round(1000 * (res["es"][1] ** 2).sum())))
