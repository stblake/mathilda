#!/usr/bin/env python3
"""Experiment 120 -- SparseArray against scipy.sparse.  See sparse_arrays.m.

Baseline: scipy.sparse CSR arrays (csr_array), the standard compressed sparse
row container; element-wise operations are scipy's own (`+`, `.multiply`,
`.power`, `.sin`, `.sqrt`), all of which preserve sparsity as SparseArray does.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness import bench, bench_once, check, require

import numpy as np
import scipy.sparse as sp

require(["scipy.sparse:csr_array", "numpy"])

n, nnz = 100000, 1000000
kk = np.arange(1, nnz + 1, dtype=np.int64)
ii = (48271 * kk) % 2147483647 % n
jj = (69621 * kk) % 2147483647 % n
jj2 = (69621 * kk + 7 * (kk % 2)) % 2147483647 % n
vv = (kk % 97) / 10.0 + 0.1


def digest(s):
    return int(round(1000 * float(s.data.sum())))


res = {}

# ---- build ------------------------------------------------------------------

bench("build: 10^6 entries, 10^5 x 10^5 (coordinates -> rows)",
      lambda: res.__setitem__("A", sp.csr_array((vv, (ii, jj)), shape=(n, n))))
A = res["A"]
check("build: 10^6 entries, 10^5 x 10^5 (coordinates -> rows)", digest(A))
B = sp.csr_array((vv, (ii, jj2)), shape=(n, n))

# ---- arithmetic -------------------------------------------------------------

for label, fn in [("arith: A + B", lambda: A + B),
                  ("arith: A B (element-wise)", lambda: A.multiply(B)),
                  ("arith: 2.5 A", lambda: 2.5 * A),
                  ("arith: A^2 (element-wise)", lambda: A.power(2)),
                  ("arith: A + A", lambda: A + A),
                  ("fn: Sin[A]", lambda: A.sin()),
                  ("fn: Sqrt[A]", lambda: A.sqrt())]:
    bench(label, lambda fn=fn: res.__setitem__("r", fn()))
    check(label, digest(sp.csr_array(res["r"])))

# ---- access -----------------------------------------------------------------

bench("access: row 500 x100", lambda: [A[[499], :] for _ in range(100)])
check("access: row 500 x100", int(round(1000 * float(A[[499], :].sum()))))
bench("access: element (500, 600) x1000", lambda: [A[499, 599] for _ in range(1000)])


def explicit():
    c = A.tocoo()
    return np.stack([c.row + 1, c.col + 1], axis=1), c.data


bench("access: explicit positions and values", lambda: res.__setitem__("e", explicit()))
check("access: explicit positions and values", int(res["e"][0].sum()))

# ---- dense <-> sparse -------------------------------------------------------

i = np.arange(1, 2001).reshape(-1, 1)
j = np.arange(1, 2001).reshape(1, -1)
dense = np.where((i + 3 * j) % 97 == 0, 1.5 + (i * j) % 7, 0.0)
bench("dense: 2000 x 2000 dense -> sparse", lambda: res.__setitem__("D", sp.csr_array(dense)))
D = res["D"]
check("dense: 2000 x 2000 dense -> sparse", digest(D))
bench("dense: 2000 x 2000 sparse -> dense", lambda: res.__setitem__("dd", D.toarray()))
check("dense: 2000 x 2000 sparse -> dense", int(round(1000 * float(res["dd"].sum()))))

# ---- products ---------------------------------------------------------------

x = (np.arange(1, n + 1) % 13) / 7.0
bench("dot: A . x (SpMV) x100", lambda: [A @ x for _ in range(100)])
check("dot: A . x (SpMV) x100", int(round(1000 * float((A @ x).sum()))))
bench("dot: x . A x10", lambda: [x @ A for _ in range(10)])
check("dot: x . A x10", int(round(1000 * float((x @ A).sum()))))
bench("dot: A . A (SpGEMM)", lambda: res.__setitem__("r", A @ A))
check("dot: A . A (SpGEMM)", digest(res["r"]))
bench("struct: Transpose[A]", lambda: res.__setitem__("At", A.T.tocsr()))
At = res["At"]
check("struct: Transpose[A]", digest(At))
bench("dot: A . Transpose[A] (SpGEMM)", lambda: res.__setitem__("r", A @ At))
check("dot: A . Transpose[A] (SpGEMM)", digest(res["r"]))

# ---- reductions -------------------------------------------------------------

import scipy.sparse.linalg as spla

for label, fn, dig in [
        ("reduce: Total[A] (column sums)", lambda: A.sum(axis=0), lambda r: int(round(1000 * float(r.sum())))),
        ("reduce: Total[A, {2}] (row sums)", lambda: A.sum(axis=1), lambda r: int(round(1000 * float(r.sum())))),
        ("reduce: Total[A, 2]", lambda: A.sum(), lambda r: int(round(1000 * float(r)))),
        ("reduce: Max[A]", lambda: A.max(), lambda r: int(round(1000 * float(r)))),
        ("reduce: Norm[A, Frobenius]", lambda: spla.norm(A, "fro"), lambda r: int(round(1000 * float(r))))]:
    bench(label, lambda fn=fn: res.__setitem__("r", fn()))
    check(label, dig(res["r"]))

# ---- assembly ---------------------------------------------------------------

bench("struct: Join[A, B] (vstack)", lambda: res.__setitem__("r", sp.vstack([A, B], format="csr")))
check("struct: Join[A, B] (vstack)", digest(res["r"]))
bench("struct: ArrayFlatten 2x2 blocks", lambda: res.__setitem__("r", sp.block_array([[A, B], [B, A]], format="csr")))
check("struct: ArrayFlatten 2x2 blocks", digest(res["r"]))
bench("struct: IdentityMatrix[10^6, SparseArray]", lambda: res.__setitem__("r", sp.identity(10**6, format="csr")))
check("struct: IdentityMatrix[10^6, SparseArray]", int(round(float(res["r"].data.sum()))))

# ---- 2-D Poisson: 5-point Laplacian on a 300 x 300 grid ---------------------

m = 300


def build_l():
    t1 = sp.diags_array([-1.0, 2.0, -1.0], offsets=[-1, 0, 1], shape=(m, m), format="csr")
    i1 = sp.identity(m, format="csr")
    return (sp.kron(t1, i1) + sp.kron(i1, t1)).tocsr()


bench("poisson: build L (Band, KroneckerProduct)", lambda: res.__setitem__("L", build_l()))
L = res["L"]
check("poisson: build L (Band, KroneckerProduct)", digest(L))
bb = (np.arange(1, m * m + 1) % 7) - 3.0
bench("poisson: LinearSolve[L, b] direct", lambda: res.__setitem__("sol", spla.spsolve(L.tocsc(), bb)))
check("poisson: LinearSolve[L, b] direct", int(round(1000 * float(res["sol"].sum()))))
bench("poisson: LinearSolve[L, b] CG 1e-10", lambda: res.__setitem__("cg", spla.cg(L, bb, rtol=1e-10)[0]))
check("poisson: LinearSolve[L, b] CG 1e-10", int(round(float(res["cg"].sum()))))
# ---- native sparse solvers (v0.6180-v0.6188) on the same L -------------------
#
# Direct methods: scipy has no CHOLMOD, so both direct rows are SuperLU --
# "Multifrontal" is spsolve (COLAMD), "Cholesky" is splu with the symmetric
# MMD_AT_PLUS_A ordering, the closest thing to a symmetric factorisation scipy
# ships.  Krylov rows use scipy's cg / bicgstab / gmres (restart 50) at the same
# 1e-10 relative tolerance.  scipy has no ILU0 or ILUTP: its only incomplete
# factorisation is SuperLU's threshold ILU (spilu), used here in natural order
# with no pivoting -- drop 1e-2 / fill 1 stands in for ILU0, drop 1e-4 / fill 5
# for ILUT and ILUTP.  That ILU is not symmetric, and CG preconditioned with it
# stalls (measured: no convergence in 1000 iterations), so the CG row runs
# scipy's cg unpreconditioned -- the honest best scipy offers for that row.

def nrm10(x):
    return int(round(10 * float(np.linalg.norm(x))))


def spilu_op(kind):
    kw = dict(drop_tol=1e-2, fill_factor=1) if kind == "ILU0" else dict(drop_tol=1e-4, fill_factor=5)
    P = spla.spilu(L.tocsc(), permc_spec="NATURAL", diag_pivot_thresh=0.0, **kw)
    return spla.LinearOperator(L.shape, P.solve)


def krylov(meth, kind):
    extra = {"restart": 50} if meth == "gmres" else {}
    M = spilu_op(kind) if kind else None
    return getattr(spla, meth)(L, bb, rtol=1e-10, M=M, maxiter=5000, **extra)[0]


Lc = L.tocsc()
bench("poisson: LinearSolve Multifrontal", lambda: res.__setitem__("sol", spla.spsolve(Lc, bb)))
check("poisson: LinearSolve Multifrontal", int(round(1000 * float(res["sol"].sum()))))
bench("poisson: LinearSolve Cholesky",
      lambda: res.__setitem__("sol", spla.splu(Lc, permc_spec="MMD_AT_PLUS_A").solve(bb)))
check("poisson: LinearSolve Cholesky", int(round(1000 * float(res["sol"].sum()))))

for label, meth, kind in [("poisson: Krylov CG + ILU0 1e-10", "cg", None),
                          ("poisson: Krylov BiCGSTAB + ILU0 1e-10", "bicgstab", "ILU0"),
                          ("poisson: Krylov BiCGSTAB + ILUT 1e-10", "bicgstab", "ILUT"),
                          ("poisson: Krylov GMRES + ILUT 1e-10", "gmres", "ILUT"),
                          ("poisson: Krylov GMRES + ILUTP 1e-10", "gmres", "ILUTP")]:
    # One timed run, no warm-up: scipy's spilu-preconditioned rows cost 8-20 s
    # each, and a warm-up would spend another minute of the file budget.
    bench_once(label, lambda meth=meth, kind=kind: res.__setitem__("x", krylov(meth, kind)))
    check(label, nrm10(res["x"]))


def lsf20():
    lu = spla.splu(Lc)
    for k in range(1, 21):
        res["x"] = lu.solve(bb + k)
    res["lu"] = lu


bench("poisson: LinearSolveFunction, factor + 20 solves", lsf20)
check("poisson: LinearSolveFunction, factor + 20 solves",
      int(round(1000 * float(res["lu"].solve(bb + 20.0).sum()))))

u0 = np.sin(np.arange(1, m * m + 1) / 100.0)
bench("poisson: MatrixExp[-0.5 L, u] (action)",
      lambda: res.__setitem__("y", spla.expm_multiply(-0.5 * L, u0)))
check("poisson: MatrixExp[-0.5 L, u] (action)", int(round(1000 * float(res["y"].sum()))))

Lz = (L + 0.5j * sp.identity(m * m, format="csr")).tocsc()
bench("poisson: complex shift LinearSolve[L + 0.5 I, b]",
      lambda: res.__setitem__("z", spla.spsolve(Lz, bb)))
check("poisson: complex shift LinearSolve[L + 0.5 I, b]",
      "{%d, %d}" % (round(1000 * res["z"].sum().real), round(1000 * res["z"].sum().imag)))


def lap2d(k):
    t = sp.diags_array([-1.0, 2.0, -1.0], offsets=[-1, 0, 1], shape=(k, k), format="csr")
    e = sp.identity(k, format="csr")
    return (sp.kron(t, e) + sp.kron(e, t)).tocsc()


def sparse_logdet(a):
    """log det of a sparse matrix from its SuperLU factors (det > 0 here)."""
    lu = spla.splu(a)
    return float(np.log(np.abs(lu.U.diagonal())).sum())


L30 = lap2d(30)
bench("sparse: Det of 900x900 Laplacian / 4", lambda: res.__setitem__("d", sparse_logdet(L30 / 4.0)))
check("sparse: Det of 900x900 Laplacian / 4", int(round(1e6 * res["d"])))
bench("sparse: Inverse of 900x900 Laplacian", lambda: res.__setitem__("inv", spla.inv(L30)))
check("sparse: Inverse of 900x900 Laplacian", int(round(1000 * float(res["inv"].sum()))))

T1k = sp.diags_array([-1.0, 2.0, -1.0], offsets=[-1, 0, 1], shape=(1000, 1000), format="csc")
bench("sparse: MatrixExp[-0.1 T] 1-D Laplacian 1000x1000",
      lambda: res.__setitem__("e", spla.expm(-0.1 * T1k)))
check("sparse: MatrixExp[-0.1 T] 1-D Laplacian 1000x1000", int(round(1000 * float(res["e"].sum()))))

bench("poisson: 6 largest eigenvalues", lambda: res.__setitem__("ev", spla.eigsh(L, k=6, which="LM", return_eigenvectors=False)))
check("poisson: 6 largest eigenvalues", int(round(1e6 * float(res["ev"].sum()))))
bench("poisson: 6 smallest eigenvalues", lambda: res.__setitem__("ev", spla.eigsh(L, k=6, sigma=0, which="LM", return_eigenvectors=False)))
check("poisson: 6 smallest eigenvalues", int(round(1e6 * float(res["ev"].sum()))))


# ---- sparse eigensystem / singular values with a count (ARPACK) --------------
# H = L + V and the graded G, not L and Dm: see the .m half for why.

_i = np.arange(1, m + 1)[:, None]
_j = np.arange(1, m + 1)[None, :]
Hpot = (L + sp.diags_array(0.001 * (_i ** 2 + 2 * _j ** 2).ravel().astype(float))).tocsr()


def es6():
    w, v = spla.eigsh(Hpot, k=6, which="LM")
    res["es"] = (w, v)


bench("schrodinger: Eigensystem[L + V, 6]", es6, reps=1)
check("schrodinger: Eigensystem[L + V, 6]",
      "{%d, %d}" % (round(1e6 * res["es"][0].sum()), round(1000 * (res["es"][1] ** 2).sum())))
g = 20000
G = sp.diags_array([100.0 / np.arange(1, g + 1), 0.1 * np.ones(g - 1), 0.2 * np.ones(g - 3)],
                   offsets=[0, 1, -3], format="csr")
bench("graded 20000: SingularValueList[G, 6]",
      lambda: res.__setitem__("sv", spla.svds(G, k=6, return_singular_vectors=False)), reps=1)
check("graded 20000: SingularValueList[G, 6]", int(round(1e6 * float(res["sv"].sum()))))
