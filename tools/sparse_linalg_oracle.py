#!/usr/bin/env python3
"""
sparse_linalg_oracle.py -- differential stress test of linear algebra on
SparseArray, against three independent oracles.

Every case is one (matrix, operation) pair.  The matrix comes from a seeded
LCG in one of fourteen classes (SPD / Hermitian positive definite, symmetric
indefinite, nonsymmetric, Hermitian, general complex, banded, block-diagonal,
singular, tall, wide, ill-conditioned, all-zero, nonzero background, 1x1) and
one of seven entry kinds (int, rational, machine real, machine complex,
30-digit MPFR real and complex, symbolic).  The operation is one of ~60
templates over LinearSolve (and its Methods), LinearSolveFunction, Det,
Inverse, MatrixRank, NullSpace, RowReduce, LeastSquares, PseudoInverse, the
eigen- and singular-value heads, the LU / QR / Cholesky / Schur / Hessenberg /
Jordan decompositions, CharacteristicPolynomial, MatrixExp / Log / Function /
Power, Minors, Adjugate, Orthogonalize, Norm, Tr, Dot, Transpose,
KroneckerProduct and the matrix predicates.

Oracles (each one optional, each batched into ONE process per run):
  reference  wolframscript, one program for every case (the licence may allow
             only one kernel at a time, so it is never run in parallel; a run
             that cannot get a licence is retried once).
  scipy      python3.11 + numpy/scipy (--python, $HPC_PYTHON), one process.
             Reads each matrix from a MatrixMarket file.
  maxima     exact rational / symbolic answers for small int, rational and
             symbolic cases (tools/data/sparse_linalg_oracle.mac), one batch.

Checks per case (one TSV row per check; `oracle` column):
  mathilda  the run itself: CRASH / TIMEOUT, or UNEVAL when Mathilda leaves the
            call unevaluated and an oracle answers it (a head still being
            built: never counted as FAIL).
  sig       the shape signature of the answer (S SparseArray, L full List,
            x scalar, U unevaluated, E $Failed, {..} ragged, Head other; the
            sparse_linalg_semantics.py convention) equals the reference's for
            the same call.
  msg       every message tag the reference emits for the call, Mathilda emits
            too.
  ref / scipy / maxima
            the answer's value.  Exact kinds compare exactly (algebraic values
            at 50 digits, to 1e-30); symbolic kinds are compared exactly at the
            rational point xa = 3/7, xb = -5/11 (a stronger test than
            Simplify[a - b] == 0, and decidable); machine kinds to
            rtol 1e-8 max(1, cond) (capped at 1e-3), MPFR kinds to the
            reference's own precision (about 25 digits).  Spectra are compared
            as multisets; non-unique results are not compared by value.
  inv       an invariant computed INSIDE Mathilda for results that are not
            unique: L.U == P.A, Q^H.R == A with Q unitary, A.V == U.S with U, V
            orthonormal, A.v == lambda v, C^H.C == A, P.H.P^H == A with H
            Hessenberg, Q.T.Q^H == A with T quasi-triangular, S.J == A.S,
            MatrixExp[MatrixLog[A]] == A, NullSpace annihilated, A.x == b for
            a consistent singular / rectangular system, ...  The reference
            computes the same invariant; when the reference fails it too, the
            row is SKIP (the invariant does not apply).
  self      Normal[f[S]] == f[Normal[S]] inside Mathilda, except where the
            reference also answers the two differently, or DIFF_OK says so.
  order     Eigenvalues / SingularValueList come in non-increasing |value|.

Verdicts: PASS FAIL SKIP UNEVAL CRASH TIMEOUT.  A case's verdict is the worst
of its rows (CRASH > TIMEOUT > FAIL > UNEVAL > PASS > SKIP).  EXEMPT turns a
matching FAIL row into SKIP with its reason.

Usage:
    python3 tools/sparse_linalg_oracle.py                   # 300 cases, seed 1
    python3 tools/sparse_linalg_oracle.py -n 1000 --seed 7 --large 8
    python3 tools/sparse_linalg_oracle.py --only 'Eigen*,Det' --kinds real
    python3 tools/sparse_linalg_oracle.py --quick           # per-push ctest
    python3 tools/sparse_linalg_oracle.py --record          # refresh --quick data
    python3 tools/sparse_linalg_oracle.py --valgrind -n 40
    python3 tools/sparse_linalg_oracle.py --repro 17        # print case 17

--quick replays a fixed case set against oracle answers recorded by --record
into tests/data/sparse_linalg_oracle/recorded.json, so it needs neither
wolframscript, scipy nor maxima.

Exit status: 1 on any FAIL or CRASH (unless --no-fail), 2 on a usage error or
a stale recording, 0 otherwise.
"""

import argparse
import fnmatch
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal, getcontext
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BINARY = os.environ.get("MATHILDA", os.path.join(ROOT, "Mathilda"))
MAXIMA_DRIVER = os.path.join(ROOT, "tools", "data", "sparse_linalg_oracle.mac")
RECORD_DIR = os.path.join(ROOT, "tests", "data", "sparse_linalg_oracle")
RECORD_FILE = os.path.join(RECORD_DIR, "recorded.json")
PYTHON = os.environ.get("HPC_PYTHON", "/usr/local/bin/python3.11")
WOLFRAMSCRIPT = os.environ.get("WOLFRAMSCRIPT",
                               "/Applications/Mathematica.app/Contents/MacOS/wolframscript")
MAXIMA = os.environ.get("MAXIMA", shutil.which("maxima") or "/usr/local/bin/maxima")

QUICK_SEED = 20261105
QUICK_N = 150
GEN_VERSION = 1     # bump when the generator changes what a seed produces

# --------------------------------------------------------------------------
# Kinds and classes
# --------------------------------------------------------------------------

KINDS = ["int", "rational", "real", "complex", "mpfr", "mpfrc", "sym"]
COMPLEX_KINDS = ("complex", "mpfrc")
GROUP = {"int": "exact", "rational": "exact", "sym": "sym", "real": "machine",
         "complex": "machine", "mpfr": "mpfr", "mpfrc": "mpfr"}
# the prelude's $SLmode for a group
MODE = {"exact": "exact", "sym": "sym", "machine": "machine", "mpfr": "mpfr"}

CLASSES = ["spd", "symind", "nonsym", "herm", "cplx", "band", "block", "sing",
           "tall", "wide", "ill", "zero", "bg", "one"]
# classes whose square matrix is nonsingular by construction (strict diagonal
# dominance, or L.D.R with unit-triangular L, R)
NONSING = {"spd", "symind", "nonsym", "herm", "cplx", "band", "block", "ill", "bg", "one"}
HPD = {"spd", "herm"}          # Hermitian positive definite (spd with a complex kind is HPD)
SYMMETRIC = {"spd", "symind", "herm"}
LARGE_CLASSES = ["spd", "band", "nonsym", "symind", "herm"]

# ill-conditioned matrices have cond ~ 1e10 by construction; used when scipy
# is not there to measure it
DEFAULT_COND = {"ill": 1e10, "sing": float("inf"), "zero": float("inf")}

# --------------------------------------------------------------------------
# Operations
# --------------------------------------------------------------------------
# Template variables (WL): A the SparseArray, n x m; x0 = (-1)^i i/(i+1),
# i = 1..m (numericised for machine / MPFR kinds); b = A.x0 (consistent for
# every class); Bm = {A.x0, A.Reverse[x0]} as columns; t a free symbol.
#
#   shape   'sq' square only, 'any' any shape
#   classes allowed classes (None = all)
#   nmax    max n per kind group; a group that is missing is not run
#   cmp     value comparison: 'val' normwise, 'det', 'set' multiset,
#           'absset' multiset of |x|, 'count'/'bool' exact, None (invariant
#           only).  When uniq(case) is false the value is not compared.
#   val     WL expression of r (the answer) that is compared
#   inv     invariant (key of INVARIANTS) evaluated inside Mathilda
#   scipy / maxima   oracle keys (None = no such oracle for the op)
#   cs      tolerance scales with the condition number
#   tol     tolerance multiplier (Krylov, characteristic polynomial)
#   large   also run in the --large tier (machine kinds, scipy only)

E, S_, P, M = "exact", "sym", "mpfr", "machine"


class Op:
    def __init__(self, oid, tpl, shape="sq", classes=None, nmax=None, cmp="val", val="r",
                 inv=None, inv_cond=False, scipy=None, maxima=None, cs=True, tol=1.0,
                 large=False, uniq=None, order=False, self_=True):
        self.id, self.tpl, self.shape = oid, tpl, shape
        self.classes = set(classes) if classes else None
        self.nmax = nmax or {E: 12, S_: 4, P: 12, M: 400}
        self.cmp, self.val, self.inv, self.inv_cond = cmp, val, inv, inv_cond
        self.scipy, self.maxima, self.cs, self.tol = scipy, maxima, cs, tol
        self.large, self.uniq, self.order, self.self_ = large, uniq, order, self_
        self.head = "Dot" if tpl.startswith("A . ") else tpl.split("[", 1)[0]

    def unique(self, case):
        return True if self.uniq is None else self.uniq(case)


def nonsingular(case):
    return case.cls in NONSING and case.rows == case.cols


def full_col_rank(case):
    return case.cls in NONSING or case.cls == "tall"


LS_N = {E: 100, S_: 4, P: 60, M: 400}
OPS = [
    # ---- LinearSolve -------------------------------------------------------
    Op("LinearSolve", "LinearSolve[A, b]", shape="any", nmax=LS_N, inv="ls",
       scipy="solve", maxima="solve", uniq=nonsingular, large=True),
    Op("LinearSolve-sparseRhs", "LinearSolve[A, SparseArray[b]]", shape="any", inv="ls",
       scipy="solve", maxima="solve", uniq=nonsingular),
    Op("LinearSolve-matRhs", "LinearSolve[A, Bm]", shape="any", inv="lsm", scipy="solveM",
       maxima="solveM", uniq=nonsingular),
    Op("LinearSolve-Multifrontal", "LinearSolve[A, b, Method -> \"Multifrontal\"]", inv="ls",
       scipy="solve", uniq=nonsingular, nmax={E: 12, P: 30, M: 400}),
    Op("LinearSolve-Cholesky", "LinearSolve[A, b, Method -> \"Cholesky\"]", inv="ls",
       classes=["spd", "herm", "symind"], scipy="solve", uniq=nonsingular,
       nmax={E: 12, P: 30, M: 400}),
    Op("LinearSolve-Krylov", "LinearSolve[A, b, Method -> \"Krylov\"]", inv="ls", tol=1e3,
       classes=NONSING - {"ill", "one"}, scipy="solve", uniq=nonsingular, large=True,
       nmax={P: 30, M: 400}),
    Op("LinearSolve-CG", "LinearSolve[A, b, Method -> {\"Krylov\", Method -> "
       "\"ConjugateGradient\", \"Preconditioner\" -> \"ILU0\"}]", inv="ls", tol=1e3,
       classes=HPD, scipy="solve", uniq=nonsingular, large=True, nmax={P: 30, M: 400}),
    Op("LinearSolve-GMRES", "LinearSolve[A, b, Method -> {\"Krylov\", Method -> \"GMRES\"}]",
       inv="ls", tol=1e3, classes=NONSING - {"ill", "one"}, scipy="solve", uniq=nonsingular,
       large=True, nmax={P: 30, M: 400}),
    Op("LinearSolve-BiCGSTAB", "LinearSolve[A, b, Method -> {\"Krylov\", Method -> "
       "\"BiCGSTAB\", \"Preconditioner\" -> \"ILU0\"}]", inv="ls", tol=1e3,
       classes=NONSING - {"ill", "one"}, scipy="solve", uniq=nonsingular, large=True,
       nmax={P: 30, M: 400}),
    Op("LinearSolveFunction", "LinearSolve[A][b]", inv="ls", scipy="solve", maxima="solve",
       uniq=nonsingular, large=True, nmax=LS_N),
    Op("LinearSolveFunction-T", "LinearSolve[A][b, \"T\"]", inv="lsT", scipy="solveT",
       maxima="solveT", uniq=nonsingular),
    # ---- determinant, inverse, rank, null space ------------------------------
    Op("Det", "Det[A]", cmp="det", scipy="det", maxima="det", nmax=LS_N),
    Op("Inverse", "Inverse[A]", scipy="inv", maxima="inv", uniq=nonsingular,
       nmax={E: 12, S_: 3, P: 30, M: 400}),
    Op("MatrixRank", "MatrixRank[A]", shape="any", cmp="count", scipy="rank", maxima="rank",
       nmax={E: 100, S_: 4, P: 30, M: 400}),
    Op("NullSpace", "NullSpace[A]", shape="any", cmp="count", val="Length[r]", inv="null",
       scipy="nullcount", maxima="nullcount", nmax={E: 12, S_: 4, P: 30, M: 200}),
    Op("RowReduce", "RowReduce[A]", shape="any", nmax={E: 12, S_: 3, P: 12, M: 100}),
    Op("LeastSquares", "LeastSquares[A, b]", shape="any", inv="lsq", scipy="lstsq",
       nmax={E: 12, S_: 3, P: 12, M: 200}),
    Op("PseudoInverse", "PseudoInverse[A]", shape="any", inv="pinv", scipy="pinv",
       nmax={E: 8, S_: 3, P: 12, M: 200}),
    # ---- spectra ----------------------------------------------------------
    Op("Eigenvalues", "Eigenvalues[A]", cmp="set", scipy="eig", maxima="eig", order=True,
       nmax={E: 4, P: 12, M: 200}),
    Op("Eigenvalues-k", "Eigenvalues[A, {k}]", cmp="absset", scipy="eigk", large=True,
       order=True, nmax={E: 4, P: 12, M: 400}),
    Op("Eigenvalues-mk", "Eigenvalues[A, -{k}]", cmp="absset", scipy="eigmk", large=True,
       nmax={E: 4, P: 12, M: 400}),
    Op("Eigenvectors", "Eigenvectors[A]", cmp="count", val="Length[r]", inv="eigvec",
       nmax={E: 3, P: 12, M: 200}),
    Op("Eigenvectors-k", "Eigenvectors[A, {k}]", cmp="count", val="Length[r]", inv="eigvec",
       large=True, nmax={E: 3, P: 12, M: 400}),
    Op("Eigensystem", "Eigensystem[A]", cmp="set", val="r[[1]]", inv="eigsys", scipy="eig",
       nmax={E: 3, P: 12, M: 200}),
    Op("Eigensystem-k", "Eigensystem[A, {k}]", cmp="absset", val="r[[1]]", inv="eigsys",
       scipy="eigk", large=True, nmax={E: 3, P: 12, M: 400}),
    Op("SingularValueList", "SingularValueList[A]", shape="any", cmp="set", scipy="svl",
       order=True, cs=False, nmax={E: 4, P: 12, M: 200}),
    Op("SingularValueList-k", "SingularValueList[A, {k}]", shape="any", cmp="set",
       scipy="svlk", cs=False, large=True, nmax={E: 4, P: 12, M: 400}),
    Op("SingularValueDecomposition", "SingularValueDecomposition[A]", shape="any", cmp="set",
       val="Diagonal[r[[2]]]", inv="svd", scipy="svdall", cs=False, nmax={E: 3, P: 12, M: 200}),
    Op("SingularValueDecomposition-k", "SingularValueDecomposition[A, {k}]", shape="any",
       cmp="set", val="Diagonal[r[[2]]]", inv="svd", scipy="svdk", cs=False,
       nmax={E: 3, P: 12, M: 200}),
    # ---- decompositions ---------------------------------------------------
    Op("LUDecomposition", "LUDecomposition[A]", cmp=None, inv="lu",
       nmax={E: 12, S_: 3, P: 30, M: 400}),
    Op("QRDecomposition", "QRDecomposition[A]", shape="any", cmp=None, inv="qr",
       nmax={E: 6, S_: 2, P: 12, M: 200}),
    Op("CholeskyDecomposition", "CholeskyDecomposition[A]", classes=["spd", "herm", "symind"],
       inv="chol", scipy="chol", uniq=lambda c: c.cls in HPD, nmax={E: 8, S_: 2, P: 30, M: 400}),
    Op("SchurDecomposition", "SchurDecomposition[A]", cmp=None, inv="schur",
       nmax={E: 4, P: 12, M: 200}),
    Op("HessenbergDecomposition", "HessenbergDecomposition[A]", cmp=None, inv="hess",
       nmax={E: 4, P: 12, M: 200}),
    Op("JordanDecomposition", "JordanDecomposition[A]", cmp=None, inv="jordan",
       classes=["block", "band", "one", "sing", "spd", "nonsym", "zero"], nmax={E: 3, S_: 2}),
    Op("CharacteristicPolynomial", "CharacteristicPolynomial[A, t]", val="CoefficientList[r, t]",
       scipy="charpoly", maxima="charpoly", tol=1e3, nmax={E: 10, S_: 3, P: 12, M: 12}),
    # ---- matrix functions -------------------------------------------------
    Op("MatrixExp", "MatrixExp[A]", scipy="expm", nmax={E: 2, P: 8, M: 200}),
    Op("MatrixExp-v", "MatrixExp[A, b]", scipy="expmv", large=True, nmax={E: 2, P: 8, M: 400}),
    Op("MatrixLog", "MatrixLog[A]", classes=HPD, inv="log", inv_cond=True, scipy="logm",
       nmax={E: 2, P: 8, M: 200}),
    Op("MatrixFunction-Sin", "MatrixFunction[Sin, A]", scipy="sinm", nmax={E: 2, P: 8, M: 100}),
    Op("MatrixPower-2", "MatrixPower[A, 2]", scipy="pow2", maxima="pow2", cs=False,
       nmax={E: 12, S_: 4, P: 30, M: 400}),
    Op("MatrixPower-m1", "MatrixPower[A, -1]", scipy="inv", maxima="inv", uniq=nonsingular,
       nmax={E: 8, S_: 3, P: 30, M: 200}),
    Op("MatrixPower-half", "MatrixPower[A, 1/2]", classes=HPD, inv="sqrt", inv_cond=True,
       scipy="sqrtm", nmax={E: 2, P: 8, M: 200}),
    # ---- minors, adjugate, orthogonalisation ------------------------------
    Op("Minors", "Minors[A, Min[2, n, m]]", shape="any", nmax={E: 6, S_: 3, P: 6, M: 8}),
    Op("Adjugate", "Adjugate[A]", scipy="adj", maxima="adj", nmax={E: 6, S_: 3, P: 8, M: 30}),
    Op("Orthogonalize", "Orthogonalize[A]", shape="any", inv="ortho",
       nmax={E: 4, S_: 2, P: 8, M: 30}),
    # ---- norms and structure ------------------------------------------------
    Op("Norm-1", "Norm[A, 1]", shape="any", scipy="norm1", maxima="norm1", cs=False,
       large=True, nmax=LS_N),
    Op("Norm-2", "Norm[A]", shape="any", scipy="norm2", cs=False,
       nmax={E: 4, S_: 2, P: 12, M: 200}),
    Op("Norm-Infinity", "Norm[A, Infinity]", shape="any", scipy="normInf", maxima="normInf",
       cs=False, large=True, nmax=LS_N),
    Op("Norm-Frobenius", "Norm[A, \"Frobenius\"]", shape="any", scipy="normF", maxima="normF",
       cs=False, large=True, nmax=LS_N),
    Op("Tr", "Tr[A]", shape="any", scipy="tr", maxima="tr", cs=False, large=True, nmax=LS_N),
    Op("Dot", "A . Transpose[A]", shape="any", scipy="dotAAt", maxima="dotAAt", cs=False,
       nmax={E: 12, S_: 4, P: 30, M: 400}),
    Op("Dot-v", "A . x0", shape="any", scipy="dotv", maxima="dotv", cs=False, large=True,
       nmax=LS_N),
    Op("Transpose", "Transpose[A]", shape="any", scipy="T", maxima="T", cs=False,
       nmax={E: 12, S_: 4, P: 12, M: 200}),
    Op("KroneckerProduct", "KroneckerProduct[A, {{1, 2}, {0, -1}}]", shape="any", scipy="kron",
       maxima="kron", cs=False, nmax={E: 12, S_: 4, P: 12, M: 30}),
    Op("PositiveDefiniteMatrixQ", "PositiveDefiniteMatrixQ[A]", cmp="bool", scipy="pdq",
       nmax={E: 12, S_: 3, P: 30, M: 400}),
    Op("SymmetricMatrixQ", "SymmetricMatrixQ[A]", shape="any", cmp="bool", scipy="symq",
       maxima="symq", nmax=LS_N),
    Op("HermitianMatrixQ", "HermitianMatrixQ[A]", shape="any", cmp="bool", scipy="hermq",
       nmax=LS_N),
]
OPS_BY_ID = {o.id: o for o in OPS}
ITERATIVE = {"LinearSolve-Krylov", "LinearSolve-CG", "LinearSolve-GMRES", "LinearSolve-BiCGSTAB"}

# Invariants: WL bodies over r (the answer), A, b, Bm, x0, n, m.  Each returns
# a nonnegative residual (machine, MPFR, or a 50-digit number for exact and
# symbolic kinds; SLnz numericises).  Helpers are in SL_PRELUDE.
INVARIANTS = {
    "ls": "SLresLS[A, r, b]",
    "lsm": "SLresLS[A, r, Bm]",
    "lsT": "SLresLS[Transpose[A], r, b]",
    "lsq": "Module[{a = SLnz[A], x = SLnz[r], y = SLnz[b]}, "
           "SLrr[ConjugateTranspose[a] . (a . x - y), SLnrm[a]^2 SLnrm[x] + SLnrm[a] SLnrm[y]]]",
    "pinv": "Module[{a = SLnz[A], x = SLnz[r]}, Max[SLrr[a . x . a - a, SLnrm[a]^2 SLnrm[x]], "
            "SLrr[x . a . x - x, SLnrm[x]^2 SLnrm[a]]]]",
    "null": "Module[{a = SLnz[A], v = SLnz[r]}, If[v === {}, 0, "
            "SLrr[a . Transpose[v], SLnrm[a] SLnrm[v]]]]",
    "eigvec": "Module[{a = SLnz[A], vs = SLnz[r]}, Max[Append[Map[SLeres1[a, #] &, vs], 0]]]",
    "eigsys": "Module[{a = SLnz[A], ls = SLnz[r[[1]]], vs = SLnz[r[[2]]]}, "
              "If[Length[ls] =!= Length[vs], 1, Max[Append[Table[SLeres2[a, ls[[i]], vs[[i]]], "
              "{i, Length[vs]}], 0]]]]",
    "svd": "Module[{u = SLnz[r[[1]]], s = SLnz[r[[2]]], v = SLnz[r[[3]]], a = SLnz[A]}, "
           "Max[SLrr[a . v - u . s, a], "
           "SLnrm[ConjugateTranspose[u] . u - IdentityMatrix[Dimensions[u][[2]]]], "
           "SLnrm[ConjugateTranspose[v] . v - IdentityMatrix[Dimensions[v][[2]]]]]]",
    "lu": "Module[{q = SLnz[r[[1]]], a = SLnz[A], k}, k = Length[q]; "
          "SLrr[(LowerTriangularize[q, -1] + IdentityMatrix[k]) . UpperTriangularize[q] "
          "- a[[r[[2]]]], a]]",
    "qr": "Module[{q = SLnz[r[[1]]], R = SLnz[r[[2]]], a = SLnz[A]}, "
          "Max[SLrr[ConjugateTranspose[q] . R - a, a], "
          "SLnrm[q . ConjugateTranspose[q] - IdentityMatrix[Length[q]]], "
          "SLrr[LowerTriangularize[R, -1], a]]]",
    "chol": "Module[{u = SLnz[r], a = SLnz[A]}, Max[SLrr[ConjugateTranspose[u] . u - a, a], "
            "SLrr[LowerTriangularize[u, -1], a]]]",
    "schur": "Module[{q = SLnz[r[[1]]], t = SLnz[r[[2]]], a = SLnz[A], k}, k = Length[t]; "
             "Max[SLrr[q . t . ConjugateTranspose[q] - a, a], "
             "SLnrm[q . ConjugateTranspose[q] - IdentityMatrix[k]], "
             "SLrr[LowerTriangularize[t, -2], a], "
             "SLrr[Table[Min[Abs[t[[i + 1, i]]], Abs[t[[i + 2, i + 1]]]], {i, k - 2}], a]]]",
    "hess": "Module[{p = SLnz[r[[1]]], h = SLnz[r[[2]]], a = SLnz[A], k}, k = Length[h]; "
            "Max[SLrr[p . h . ConjugateTranspose[p] - a, a], "
            "SLnrm[p . ConjugateTranspose[p] - IdentityMatrix[k]], "
            "SLrr[LowerTriangularize[h, -2], a]]]",
    "jordan": "Module[{s = SLnz[r[[1]]], j = SLnz[r[[2]]], a = SLnz[A]}, "
              "Max[SLrr[s . j - a . s, SLnrm[a] SLnrm[s]], SLrr[LowerTriangularize[j, -1], a], "
              "SLrr[UpperTriangularize[j, 2], a], "
              "If[Abs[Det[s]] > 10^-20 SLnrm[s]^Length[s], 0, 1]]]",
    "log": "Module[{x = SLnz[r], a = SLnz[A]}, SLrr[MatrixExp[x] - a, a]]",
    "sqrt": "Module[{x = SLnz[r], a = SLnz[A]}, SLrr[x . x - a, a]]",
    "ortho": "Module[{o = Select[Normal[SLnz[r]], SLnrm[#] > 10^-10 &]}, If[o === {}, 0, "
             "SLnrm[o . ConjugateTranspose[o] - IdentityMatrix[Length[o]]]]]",
}

# --------------------------------------------------------------------------
# Deliberate deviations
# --------------------------------------------------------------------------
# (op glob, class glob, kind glob, oracle glob) -> reason.  A FAIL row that
# matches becomes SKIP "EXEMPT: reason".  Every line is a semantics decision.
#
# DIFF_OK: where the REFERENCE itself answers f[S] and f[Normal[S]]
# differently, so the `self` check (Normal[f[S]] == f[Normal[S]]) does not
# apply.
DIFF_OK = {
    ("LinearSolve*", "sing", "real", "self"):
        "machine sparse singular LinearSolve: the reference leaves it unevaluated "
        "(LinearSolve::sing), dense gives a particular solution",
    ("LinearSolve*", "sing", "complex", "self"):
        "machine sparse singular LinearSolve: LinearSolve::sing vs a dense particular solution",
    ("LinearSolve*", "tall", "real", "self"):
        "sparse rectangular machine LinearSolve gives LinearSolve::sqmat; dense solves",
    ("LinearSolve*", "tall", "complex", "self"): "as above (::sqmat)",
    ("LinearSolve*", "wide", "real", "self"): "as above (::sqmat)",
    ("LinearSolve*", "wide", "complex", "self"): "as above (::sqmat)",
}
EXEMPT = dict(DIFF_OK)
EXEMPT.update({
    ("Inverse", "*", "sym", "sig"):
        "plan (Result heads): Mathilda evaluates a symbolic sparse Inverse the "
        "reference leaves unevaluated",
    ("SingularValueList*", "zero", "*", "scipy"):
        "scipy's svds/svdvals drop zero singular values of an all-zero matrix; "
        "the reference (and Mathilda) report them as 0 (checked by the ref row)",
    ("LinearSolve-CG", "herm", "*", "msg"):
        "the reference's CG with ILU0 on a complex Hermitian matrix stalls "
        "(::krymit) -- ILU0 of a Hermitian matrix is not Hermitian; Mathilda's "
        "converges to the residual tolerance (checked by the value rows)",
})


def exempt_reason(op, cls, kind, oracle):
    for (o, c, k, w), why in EXEMPT.items():
        if (fnmatch.fnmatch(op, o) and fnmatch.fnmatch(cls, c) and fnmatch.fnmatch(kind, k)
                and fnmatch.fnmatch(oracle, w)):
            return why
    return None


# --------------------------------------------------------------------------
# Reproducible LCG (independent of Python's random module across versions)
# --------------------------------------------------------------------------

class LCG:
    MASK = (1 << 64) - 1

    def __init__(self, seed):
        self.s = (seed * 0x9E3779B97F4A7C15 + 0x632BE59BD9B4E019) & self.MASK
        for _ in range(3):
            self.next()

    def next(self):
        self.s = (self.s * 6364136223846793005 + 1442695040888963407) & self.MASK
        return self.s >> 11

    def random(self):
        return self.next() / float(1 << 53)

    def randint(self, a, b):
        return a + self.next() % (b - a + 1)

    def choice(self, seq):
        return seq[self.randint(0, len(seq) - 1)]

    def shuffle(self, x):
        for i in range(len(x) - 1, 0, -1):
            j = self.randint(0, i)
            x[i], x[j] = x[j], x[i]

    def sample(self, n, k):
        """k distinct integers from range(n) (sparse Fisher-Yates)."""
        k = min(k, n)
        swap, out = {}, []
        for i in range(k):
            j = self.randint(i, n - 1)
            vi, vj = swap.get(i, i), swap.get(j, j)
            swap[j] = vi
            out.append(vj)
        return out


# --------------------------------------------------------------------------
# Matrix generation (entries exact: complex rationals as (re, im) Fractions,
# or ('sym', wl_text, maxima_text) for a symbolic entry)
# --------------------------------------------------------------------------

F = Fraction
ZERO = (F(0), F(0))


def nzint(rng, k):
    v = 0
    while v == 0:
        v = rng.randint(-k, k)
    return v


def cval(rng, kind):
    """A random nonzero entry value of the given kind."""
    def one():
        if kind in ("int", "sym"):
            return F(nzint(rng, 9))
        if kind in ("rational", "mpfr", "mpfrc"):
            return F(nzint(rng, 9), rng.randint(1, 7))
        return F(nzint(rng, 9999), 1000)
    return (one(), one() if kind in COMPLEX_KINDS else F(0))


def cadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def cscale(c, a):
    return (c * a[0], c * a[1])


def conj(a):
    return (a[0], -a[1])


def mag(a):
    return abs(a[0]) + abs(a[1])     # >= the modulus: a safe dominance bound


class Mat:
    def __init__(self, rows, cols, entries, bg=ZERO):
        self.rows, self.cols, self.E, self.bg = rows, cols, entries, bg

    def get(self, i, j):
        return self.E.get((i, j), self.bg)


def _dominate(rng, E, n, m, signs, bg=ZERO):
    """Put a diagonal on rows 0..min(n,m)-1 that strictly dominates the row."""
    rowsum, count = {}, {}
    for (i, j), v in E.items():
        if i != j:
            rowsum[i] = rowsum.get(i, F(0)) + mag(v)
            count[i] = count.get(i, 0) + 1
    for i in range(min(n, m)):
        s = rowsum.get(i, F(0))
        if bg != ZERO:     # off-diagonal positions with no rule hold the background
            s += (m - 1 - count.get(i, 0)) * mag(bg)
        d = s + rng.randint(1, 9)
        E[(i, i)] = (signs[i] * d, F(0))


def _signs(rng, n, mixed):
    if not mixed:
        return [1] * n
    s = [rng.choice([1, -1]) for _ in range(n)]
    if n >= 2 and len(set(s)) == 1:
        s[rng.randint(0, n - 1)] *= -1
    return s


def _per_row(rng, n, tier):
    if tier == "small":
        return rng.randint(0, max(1, n // 2))
    return rng.randint(1, 5)


def gen_matrix(rng, cls, kind, n, m, tier):
    """An exact matrix of class `cls`; n x m (square classes ignore m)."""
    E = {}
    bg = ZERO
    if cls in ("spd", "herm", "symind"):
        for i in range(n):
            for j in rng.sample(n, _per_row(rng, n, tier)):
                if j > i:
                    v = cval(rng, kind)
                    E[(i, j)] = v
                    E[(j, i)] = conj(v)
        _dominate(rng, E, n, n, _signs(rng, n, cls == "symind"))
        m = n
    elif cls in ("nonsym", "cplx", "bg", "sing", "tall"):
        if cls == "bg":
            bg = cval(rng, kind)
        for i in range(n):
            for j in rng.sample(n, _per_row(rng, n, tier)):
                if j != i:
                    E[(i, j)] = cval(rng, kind)
        if cls == "bg" and n >= 2:
            for _ in range(rng.randint(1, 2)):     # explicit zeros under a nonzero background
                i, j = rng.randint(0, n - 1), rng.randint(0, n - 1)
                if i != j:
                    E[(i, j)] = ZERO
        _dominate(rng, E, n, n, _signs(rng, n, True), bg)
        m = n
        if cls == "sing":
            E = _make_singular(rng, E, n, kind)
        elif cls == "tall":
            extra = rng.randint(1, max(1, n // 2 + 1))
            for i in range(n, n + extra):
                for j in rng.sample(n, max(1, _per_row(rng, n, tier))):
                    E[(i, j)] = cval(rng, kind)
            perm = list(range(n + extra))
            rng.shuffle(perm)
            E = {(perm[i], j): v for (i, j), v in E.items()}
            n, m = n + extra, n
    elif cls == "wide":
        T = gen_matrix(rng, "tall", kind, m, m, tier)
        return Mat(T.cols, T.rows, {(j, i): v for (i, j), v in T.E.items()})
    elif cls == "band":
        w = rng.choice([1, 2])
        for i in range(n):
            for j in range(max(0, i - w), min(n, i + w + 1)):
                if j != i and rng.random() < 0.85:
                    E[(i, j)] = cval(rng, kind)
        _dominate(rng, E, n, n, _signs(rng, n, rng.random() < 0.5))
        m = n
    elif cls == "block":
        i0 = 0
        while i0 < n:
            sz = min(n - i0, rng.randint(1, 4))
            for i in range(i0, i0 + sz):
                for j in range(i0, i0 + sz):
                    if i != j and rng.random() < 0.8:
                        E[(i, j)] = cval(rng, kind)
            i0 += sz
        _dominate(rng, E, n, n, _signs(rng, n, True))
        m = n
    elif cls == "ill":
        # A = L.D.R, L unit lower bidiagonal, R unit upper bidiagonal,
        # D = 10^(-10 i/(n-1)): cond ~ 1e10.  Exact integers for int kinds.
        scale = F(10) ** 10 if kind in ("int", "sym") else F(1)
        D = [scale / F(10) ** round(10 * i / (n - 1)) for i in range(n)]
        def lr():
            # |entry| <= 1: the bidiagonal factors stay well conditioned and
            # cond(A) ~ cond(D) = 1e10 rather than compounding with n
            if rng.random() >= 0.7:
                return ZERO
            if kind in ("int", "sym"):
                return (F(rng.choice([1, -1])), F(0))
            im = F(rng.randint(-9, 9), 10) if kind in COMPLEX_KINDS else F(0)
            return (F(nzint(rng, 9), 10), im)
        L = [lr() for _ in range(n)]
        R = [lr() for _ in range(n)]
        # (L.D.R)_{ij}: L_{i,i-1} = L[i], R_{j-1,j} = R[j]
        for i in range(n):
            for j in range(max(0, i - 1), min(n, i + 2)):
                v = ZERO
                for k in (i - 1, i):
                    if k < 0 or k > j or k < j - 1:
                        continue
                    lik = (F(1), F(0)) if k == i else L[i]
                    rkj = (F(1), F(0)) if k == j else R[j]
                    v = cadd(v, cmul(cmul(lik, (D[k], F(0))), rkj))
                if v != ZERO:
                    E[(i, j)] = v
        m = n
    elif cls == "zero":
        pass
    elif cls == "one":
        E[(0, 0)] = cval(rng, kind)
        n = m = 1
    else:
        raise ValueError(cls)
    if kind == "sym" and E:
        _symbolise(rng, E, cls in SYMMETRIC)
    return Mat(n, m, E, bg)


def _make_singular(rng, E, n, kind):
    rows = {}
    for (i, j), v in E.items():
        rows.setdefault(i, {})[j] = v
    mode = rng.choice(["dup", "zero", "sum", "lowrank"] if n >= 3 else ["dup", "zero"])
    i, j = rng.sample(n, 2)
    if mode == "dup":
        c = F(rng.choice([1, -1, 2]))
        rows[j] = {col: cscale(c, v) for col, v in rows.get(i, {}).items()}
    elif mode == "zero":
        rows[j] = {}
    elif mode == "sum":
        k = [x for x in range(n) if x not in (i, j)][0]
        s = dict(rows.get(i, {}))
        for col, v in rows.get(k, {}).items():
            s[col] = cadd(s.get(col, ZERO), v)
        rows[j] = {c: v for c, v in s.items() if v != ZERO}
    else:   # rank <= 2: u v^T + w z^T over sparse vectors
        rows = {}
        for _ in range(2):
            u = {r: cval(rng, kind) for r in rng.sample(n, rng.randint(1, n))}
            v = {c: cval(rng, kind) for c in rng.sample(n, rng.randint(1, n))}
            for r, ur in u.items():
                for c, vc in v.items():
                    rows.setdefault(r, {})
                    rows[r][c] = cadd(rows[r].get(c, ZERO), cmul(ur, vc))
    return {(r, c): v for r, d in rows.items() for c, v in d.items() if v != ZERO}


SYMBOLS = [("xa", "xa"), ("xb", "xb"), ("xa + 1", "xa+1"), ("2 xb", "2*xb"), ("xa - xb", "xa-xb")]


def _symbolise(rng, E, symmetric):
    keys = sorted(E)
    for _ in range(rng.randint(1, 2)):
        i, j = rng.choice(keys)
        wl, mx = rng.choice(SYMBOLS)
        E[(i, j)] = ("sym", wl, mx)
        if symmetric:
            E[(j, i)] = ("sym", wl, mx)


# ---- rendering ------------------------------------------------------------

def frs(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def fls(x):
    s = repr(float(x))
    if "e" in s:
        mant, ex = s.split("e")
        if "." not in mant:
            mant += ".0"
        s = mant + "*^" + str(int(ex))
    return s


def wl_value(v, kind):
    if v[0] == "sym":
        return "(" + v[1] + ")"
    re_, im = v
    if kind in ("int", "rational", "sym"):
        return frs(re_)
    if kind in ("real", "complex"):
        return fls(re_) if im == 0 else f"Complex[{fls(re_)}, {fls(im)}]"
    if im == 0:
        return f"N[{frs(re_)}, 30]"
    return f"N[{frs(re_)} + ({frs(im)}) I, 30]"


def wl_matrix(mat, kind):
    rules = ", ".join(f"{{{i + 1}, {j + 1}}} -> {wl_value(v, kind)}"
                      for (i, j), v in sorted(mat.E.items()))
    tail = "" if mat.bg == ZERO else ", " + wl_value(mat.bg, kind)
    return f"SparseArray[{{{rules}}}, {{{mat.rows}, {mat.cols}}}{tail}]"


def maxima_matrix(mat):
    def mv(v):
        if v[0] == "sym":
            return "(" + v[2] + ")"
        return "(" + frs(v[0]) + ")" if v[0].denominator != 1 or v[0] < 0 else frs(v[0])
    rows = ["[" + ",".join(mv(mat.get(i, j)) for j in range(mat.cols)) + "]"
            for i in range(mat.rows)]
    return "matrix(" + ",".join(rows) + ")"


def write_mtx(path, mat, kind):
    """MatrixMarket coordinate file (real or complex), the scipy side's input."""
    cx = kind in COMPLEX_KINDS
    lines = []
    if mat.bg == ZERO:
        items = sorted(mat.E.items())
    else:
        items = [((i, j), mat.get(i, j)) for i in range(mat.rows) for j in range(mat.cols)]
    for (i, j), v in items:
        if v[0] == "sym":
            raise ValueError("symbolic entry in MTX")
        if cx:
            lines.append(f"{i + 1} {j + 1} {float(v[0])!r} {float(v[1])!r}")
        else:
            lines.append(f"{i + 1} {j + 1} {float(v[0])!r}")
    with open(path, "w") as fh:
        fh.write("%%MatrixMarket matrix coordinate " + ("complex" if cx else "real") + " general\n")
        fh.write(f"{mat.rows} {mat.cols} {len(lines)}\n")
        fh.write("\n".join(lines) + ("\n" if lines else ""))


def hadamard(mat):
    """prod of row 2-norms (|Det| bound) as a float."""
    p = 1.0
    for i in range(mat.rows):
        s = 0.0
        for j in range(mat.cols):
            v = mat.get(i, j)
            if v[0] == "sym":
                return float("inf")
            s += float(v[0]) ** 2 + float(v[1]) ** 2
        p *= math.sqrt(s)
    return p


# --------------------------------------------------------------------------
# Cases
# --------------------------------------------------------------------------

class Case:
    def __init__(self, cid, op, cls, kind, tier, mat, k):
        self.id, self.op, self.cls, self.kind, self.tier, self.mat = cid, op, cls, kind, tier, mat
        self.rows, self.cols = mat.rows, mat.cols
        self.group = GROUP[kind]
        self.tpl = op.tpl.replace("{k}", str(k))
        self.k = k
        self.wl = wl_matrix(mat, kind)
        self.fp = hashlib.sha1(f"{GEN_VERSION}|{op.id}|{self.tpl}|{cls}|{kind}|{self.wl}"
                               .encode()).hexdigest()[:16]
        self.cond = None

    @property
    def size(self):
        return f"{self.rows}x{self.cols}"


def class_ok(op, cls, kind, tier):
    if op.classes is not None and cls not in op.classes:
        return False
    if cls in ("herm", "cplx") and kind not in COMPLEX_KINDS:
        return False
    if kind == "sym" and cls in ("ill", "herm", "cplx"):
        return False
    if cls in ("tall", "wide") and op.shape != "any":
        return False
    if tier == "large":
        return cls in LARGE_CLASSES and kind in ("real", "complex")
    return True


def pick_size(rng, op, cls, kind, tier):
    nmax = op.nmax.get(GROUP[kind])
    if nmax is None:
        return None
    if cls == "one":
        return (1, 1) if tier == "small" else None
    lo, hi = {"small": (2, 12), "medium": (50, 400), "large": (2000, 20000)}[tier]
    if cls == "ill":
        lo = max(lo, 3)
    if tier != "large":
        hi = min(hi, nmax)
    if hi < lo:
        return None
    n = rng.randint(lo, hi)
    if cls == "zero" and op.shape == "any" and rng.random() < 0.4:
        return (n, rng.randint(max(1, lo - 1), hi))
    return (n, n)


def gen_cases(n, seed, medium_frac=0.2, large=0, ops=None, classes=None, kinds=None):
    master = LCG(seed)
    ops = ops or OPS
    classes = classes or CLASSES
    kinds = kinds or KINDS
    cases = []
    cycles = {}     # one op cycle per tier: the large tier draws only large-capable ops

    def next_op(pool, tier):
        cycle = cycles.setdefault(tier, [])
        if not cycle:
            cycle.extend(pool)
            master.shuffle(cycle)
        return cycle.pop()

    total = n + large
    for cid in range(total):
        tier = "large" if cid >= n else ("medium" if master.random() < medium_frac else "small")
        pool = [o for o in ops if o.large] if tier == "large" else ops
        if not pool:
            continue
        op = next_op(pool, tier)
        for _attempt in range(400):
            if _attempt and _attempt % 100 == 0:
                op = next_op(pool, tier)    # this op cannot be placed in this tier
            cls = master.choice(LARGE_CLASSES if tier == "large" else classes)
            kind = master.choice(kinds)
            if kind not in kinds or not class_ok(op, cls, kind, tier):
                continue
            dims = pick_size(master, op, cls, kind, tier)
            if dims is None:
                continue
            crng = LCG(seed * 1000003 + cid)
            mat = gen_matrix(crng, cls, kind, dims[0], dims[1], tier)
            if op.shape == "sq" and mat.rows != mat.cols:
                continue
            k = min(2, max(1, min(mat.rows, mat.cols) - (1 if tier == "large" else 0)))
            cases.append(Case(cid, op, cls, kind, tier, mat, k))
            break
    return cases


# --------------------------------------------------------------------------
# The WL program (identical for Mathilda and the reference)
# --------------------------------------------------------------------------

SL_PRELUDE = r"""
$Lsub = {xa -> 3/7, xb -> -5/11};
SLjs[s_String] := "\"" <> StringReplace[s, {"\\" -> "\\\\", "\"" -> "\\\"", "\n" -> " ", "\t" -> " "}] <> "\"";
(* one function with Which, not a family of DownValues: the encoder must not
   depend on how a system orders a general rule against specific ones *)
SLenc[x_] := Which[
  ListQ[x], "[" <> StringRiffle[Map[SLenc, x], ","] <> "]",
  x === True, "true", x === False, "false",
  IntegerQ[x], SLjs["i:" <> ToString[x]],
  Head[x] === Rational, SLjs["q:" <> ToString[Numerator[x]] <> "/" <> ToString[Denominator[x]]],
  Head[x] === Real, SLjs[If[MachineNumberQ[x], "r:", "R:"] <> ToString[x, InputForm]],
  Head[x] === Complex, "{\"c\":[" <> SLenc[Re[x]] <> "," <> SLenc[Im[x]] <> "]}",
  True, SLjs["s:" <> ToString[x, InputForm]]];
SLexactQ[y_] := And @@ Map[(IntegerQ[#] || Head[#] === Rational || # === True || # === False ||
    (Head[#] === Complex && FreeQ[#, _Real])) &, Flatten[{y}]];
(* the compared form of an answer: Normal, symbols at the rational point, and
   an exact answer that is not rational (Sqrt, Root, Exp) at 50 digits *)
SLcanon[x_] := Module[{y = x /. s_SparseArray :> Normal[s]},
  If[$SLmode === "sym", y = y /. $Lsub];
  If[($SLmode === "exact" || $SLmode === "sym") && ! SLexactQ[y], y = N[y, 50]];
  y];
(* numericised operand for an invariant: exact and symbolic at 50 digits *)
SLnz[x_] := Module[{y = x},
  If[$SLmode === "exact" || $SLmode === "sym",
    y = y /. s_SparseArray :> Normal[s];
    If[$SLmode === "sym", y = y /. $Lsub];
    y = N[y, 50]];
  y];
SLnrm[x_] := Max[Append[Abs[If[Head[x] === SparseArray, ArrayRules[x][[All, 2]], Flatten[{x}]]], 0]];
SLrr[res_, sc_] := Module[{a = SLnrm[res], c = SLnrm[sc]}, If[c == 0, a, a/c]];
SLresLS[A_, r_, b_] := Module[{a = SLnz[A], x = SLnz[r], y = SLnz[b]},
  SLrr[a . x - y, SLnrm[a] SLnrm[x] + SLnrm[y]]];
SLeres1[a_, v_] := Module[{vv = Conjugate[v] . v, l},
  If[vv == 0, 0, l = (Conjugate[v] . (a . v))/vv; SLrr[a . v - l v, (SLnrm[a] + Abs[l]) SLnrm[v]]]];
SLeres2[a_, l_, v_] := If[Conjugate[v] . v == 0, 0, SLrr[a . v - l v, (SLnrm[a] + Abs[l]) SLnrm[v]]];
SLsig[r_, hd_] := Which[
  r === $TIMEOUT, "T",
  Head[r] === hd || Head[Head[r]] === hd, "U",
  Head[r] === SparseArray, "S",
  r === $Failed || r === $Aborted, "E",
  ListQ[r] && ArrayQ[r] && FreeQ[r, _SparseArray], "L",
  ListQ[r], "{" <> StringRiffle[Map[SLsig[#, hd] &, r], ","] <> "}",
  AtomQ[r] || NumericQ[r] || MemberQ[{Plus, Times, Power, Rational, Complex}, Head[r]], "x",
  True, ToString[Head[r]]];
SLbad = {"U", "T", "E"};
SLrun[id_, mode_, hd_, F_, V_, Inv_, selfQ_, mat_] := Block[{$SLmode = mode},
 Module[{A = mat, n, m, x0, b, Bm, r, msgs, s, val = "null", inv = "null", d, dsig = "null",
     dval = "null", ds},
  n = Dimensions[A][[1]]; m = Dimensions[A][[2]];
  x0 = Table[(-1)^i i/(i + 1), {i, m}];
  If[mode === "machine", x0 = N[x0]];
  If[mode === "mpfr", x0 = N[x0, 30]];
  b = Normal[A . x0]; Bm = Transpose[{b, Normal[A . Reverse[x0]]}];
  Block[{$MessageList = {}},
    r = TimeConstrained[F[A, b, Bm, x0, n, m], $SLtl, $TIMEOUT];
    msgs = Map[ToString, $MessageList]];
  s = TimeConstrained[SLsig[r, hd], $SLtl, "T"];
  If[! MemberQ[SLbad, s],
    val = TimeConstrained[SLenc[SLcanon[V[r, A, b, Bm, x0, n, m]]], $SLtl, SLjs["s:$TIMEOUT"]];
    If[Inv =!= None,
      inv = TimeConstrained[SLenc[Inv[r, A, b, Bm, x0, n, m]], $SLtl, SLjs["s:$TIMEOUT"]]]];
  If[selfQ,
    d = TimeConstrained[F[Normal[A], b, Bm, x0, n, m], $SLtl, $TIMEOUT];
    ds = TimeConstrained[SLsig[d, hd], $SLtl, "T"];
    dsig = SLjs[ds];
    If[! MemberQ[SLbad, ds],
      dval = TimeConstrained[SLenc[SLcanon[V[d, Normal[A], b, Bm, x0, n, m]]], $SLtl, "null"]]];
  WriteString["stdout", "@@R {\"id\":", ToString[id], ",\"sig\":", SLjs[s], ",\"msgs\":[",
    StringRiffle[Map[SLjs, msgs], ","], "],\"val\":", val, ",\"inv\":", inv, ",\"dsig\":", dsig,
    ",\"dval\":", dval, "}\n"]]];
"""


def case_line(c):
    op = c.op
    F_ = f"Function[{{A, b, Bm, x0, n, m}}, {c.tpl}]"
    V = f"Function[{{r, A, b, Bm, x0, n, m}}, {op.val}]"
    I_ = f"Function[{{r, A, b, Bm, x0, n, m}}, {INVARIANTS[op.inv]}]" if op.inv else "None"
    selfq = "True" if (op.self_ and op.cmp and c.tier != "large") else "False"
    return (f'SLrun[{c.id}, "{MODE[c.group]}", {op.head}, {F_}, {V}, {I_}, {selfq}, '
            f'{c.wl}];')


def program(cases, timelimit):
    return SL_PRELUDE + f"$SLtl = {timelimit};\n" + "\n".join(case_line(c) for c in cases) + "\n"


def parse_result_lines(text, tag="@@R "):
    out = {}
    for line in text.splitlines():
        if line.startswith(tag):
            try:
                d = json.loads(line[len(tag):])
            except ValueError:
                continue
            out[d["id"]] = d
    return out


# --------------------------------------------------------------------------
# Running Mathilda (chunks, crash bisection, valgrind)
# --------------------------------------------------------------------------

MEMORY_ERROR_KINDS = ("Invalid read", "Invalid write", "Invalid free", "Mismatched free",
                      "Conditional jump", "Use of uninitialised value",
                      "Source and destination overlap")


def _run(cmd, timeout, cwd=ROOT):
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, timeout=timeout)
        return p.returncode, p.stdout.decode(errors="replace"), p.stderr.decode(errors="replace")
    except subprocess.TimeoutExpired as e:
        out = e.stdout or b""
        err = e.stderr or b""
        if isinstance(out, bytes):
            out = out.decode(errors="replace")
        if isinstance(err, bytes):
            err = err.decode(errors="replace")
        return -999, out, err + "\n[chunk timeout]"


def memory_errors(err):
    for blk in re.split(r"\n==\d+== \n", err):
        if any(k in blk for k in MEMORY_ERROR_KINDS):
            return blk[-1500:]
    return ""


def definitely_lost(err):
    m = re.search(r"definitely lost: ([0-9,]+) bytes", err)
    return int(m.group(1).replace(",", "")) if m else 0


class MathildaRunner:
    def __init__(self, binary, workdir, timelimit, valgrind=False, keep=None):
        self.binary, self.workdir, self.tl = binary, workdir, timelimit
        self.valgrind, self.keep = valgrind, keep
        self._n = 0
        self._lock = threading.Lock()
        self._leak_base = None

    def _cmd(self, path):
        cmd = [self.binary, "-file", path]
        if self.valgrind:
            cmd = ["valgrind", "--leak-check=full", "--errors-for-leak-kinds=none"] + cmd
        return cmd

    def leak_baseline(self):
        if self._leak_base is None:
            path = os.path.join(self.workdir, "leakbase.m")
            with open(path, "w") as fh:
                fh.write(SL_PRELUDE)
            _, _, err = _run(self._cmd(path), 900)
            self._leak_base = definitely_lost(err)
        return self._leak_base

    def run_chunk(self, cases):
        with self._lock:
            self._n += 1
            k = self._n
        path = os.path.join(self.workdir, f"mathilda_{k:04d}.m")
        with open(path, "w") as fh:
            fh.write(program(cases, self.tl))
        timeout = min(1800, 60 + 4 * self.tl * len(cases)) * (20 if self.valgrind else 1)
        rc, out, err = _run(self._cmd(path), timeout)
        extra = {}
        if self.valgrind:
            merr = memory_errors(err)
            if merr:
                rc, err = 99, merr
            lost = definitely_lost(err)
            base = self.leak_baseline()
            if lost > base:
                extra[cases[0].id] = f"definitely lost {lost} bytes (startup baseline {base})"
        if self.keep:
            with open(path[:-2] + ".out", "w") as fh:
                fh.write(out)
            with open(path[:-2] + ".err", "w") as fh:
                fh.write(err[-200000:])
        else:
            os.unlink(path)
        return rc, parse_result_lines(out), err, extra

    def run_all(self, cases, chunk, jobs):
        """id -> result dict; crashes bisected down to the case."""
        results, crashes, leaks = {}, {}, {}
        pending = [cases[i:i + chunk] for i in range(0, len(cases), chunk)]
        lock = threading.Lock()

        def work(part):
            rc, res, err, extra = self.run_chunk(part)
            with lock:
                results.update(res)
                leaks.update(extra)
            if rc == 0:
                return []
            rest = [c for c in part if c.id not in res]
            if not rest:
                return []
            if len(part) == 1:
                kind = "TIMEOUT" if rc == -999 else "CRASH"
                with lock:
                    crashes[rest[0].id] = (kind, f"rc={rc} " + err.strip()[-300:].replace("\n", " | "))
                return []
            if len(rest) == 1:
                return [rest]
            mid = len(rest) // 2
            return [rest[:mid], rest[mid:]]

        with ThreadPoolExecutor(max_workers=max(1, jobs)) as ex:
            while pending:
                futs = [ex.submit(work, p) for p in pending]
                pending = []
                for f in futs:
                    pending.extend(f.result())
        return results, crashes, leaks


# --------------------------------------------------------------------------
# The reference (wolframscript)
# --------------------------------------------------------------------------

def run_reference(cases, workdir, timelimit, wolframscript, log):
    if not cases:
        return {}
    path = os.path.join(workdir, "reference.m")
    with open(path, "w") as fh:
        fh.write(program(cases, timelimit))
    timeout = 300 + len(cases) * (2 * timelimit + 5)
    for attempt in (1, 2):
        t0 = time.time()
        rc, out, err = _run([wolframscript, "-file", path], timeout)
        with open(os.path.join(workdir, "reference.out"), "w") as fh:
            fh.write(out + "\n(* stderr *)\n" + err[-200000:])
        res = parse_result_lines(out)
        log(f"reference: {len(res)}/{len(cases)} answers in {time.time() - t0:.0f} s"
            + ("" if rc == 0 else f" (rc={rc})"))
        if res or attempt == 2:
            return res
        log("reference: no answers (licence busy?); retrying once in 30 s: "
            + (out + err).strip()[-200:].replace("\n", " "))
        time.sleep(30)
    return {}


# --------------------------------------------------------------------------
# scipy (one python3.11 process; this file is the worker, --scipy-worker)
# --------------------------------------------------------------------------

def scipy_jobs(cases, workdir):
    jobs = []
    mtxdir = os.path.join(workdir, "mtx")
    os.makedirs(mtxdir, exist_ok=True)
    for c in cases:
        if c.kind == "sym":
            continue
        path = os.path.join(mtxdir, f"case_{c.id}.mtx")
        write_mtx(path, c.mat, c.kind)
        key = c.op.scipy if (c.op.scipy and c.group in ("exact", "machine")
                             and c.op.unique(c)) else None
        if c.op.id == "PositiveDefiniteMatrixQ" and c.cls not in SYMMETRIC:
            key = None  # the reference's definition for a nonsymmetric matrix is its own
        jobs.append({"id": c.id, "op": key, "mtx": path, "rows": c.rows, "cols": c.cols,
                     "k": c.k, "cls": c.cls, "large": c.tier == "large",
                     "sym": c.cls in SYMMETRIC})
    return jobs


def run_scipy(cases, workdir, python, log):
    jobs = scipy_jobs(cases, workdir)
    if not jobs:
        return {}
    path = os.path.join(workdir, "scipy_jobs.json")
    with open(path, "w") as fh:
        json.dump(jobs, fh)
    t0 = time.time()
    rc, out, err = _run([python, os.path.abspath(__file__), "--scipy-worker", path],
                        600 + 30 * len(jobs))
    with open(os.path.join(workdir, "scipy.out"), "w") as fh:
        fh.write(out + err)
    res = parse_result_lines(out, "@@S ")
    log(f"scipy: {len(res)}/{len(jobs)} answers in {time.time() - t0:.0f} s"
        + ("" if rc == 0 else f" (rc={rc}: {err.strip()[-300:]})"))
    return res


def scipy_worker(path):
    """Runs under python3.11 with numpy/scipy: one @@S line per job."""
    import warnings
    import numpy as np
    import scipy.io
    import scipy.linalg as sla
    import scipy.sparse as sps
    import scipy.sparse.linalg as spla
    warnings.simplefilter("ignore")

    def enc(x):
        if isinstance(x, (bool, np.bool_)):
            return bool(x)
        if isinstance(x, (int, np.integer)):
            return f"i:{int(x)}"
        if isinstance(x, np.ndarray):
            return [enc(v) for v in x]
        if isinstance(x, (list, tuple)):
            return [enc(v) for v in x]
        if isinstance(x, (complex, np.complexfloating)):
            return {"c": [f"r:{float(x.real)!r}", f"r:{float(x.imag)!r}"]}
        return f"r:{float(x)!r}"

    with open(path) as fh:
        jobs = json.load(fh)
    K = np.array([[1.0, 2.0], [0.0, -1.0]])
    for j in jobs:
        out = {"id": j["id"]}
        try:
            A = sps.csr_matrix(scipy.io.mmread(j["mtx"]))
            m = j["cols"]
            x0 = np.array([(-1) ** i * i / (i + 1) for i in range(1, m + 1)])
            large = j["large"]
            D = None if large else A.toarray()
            if not large:
                try:
                    out["cond"] = float(np.linalg.cond(D)) if D.size else float("inf")
                except Exception:
                    out["cond"] = float("inf")
            b = A @ x0
            op, k = j["op"], j["k"]
            v = None
            if op is None:
                pass
            elif op == "solve":
                v = spla.spsolve(A.tocsc(), b) if large else np.linalg.solve(D, b)
            elif op == "solveM":
                v = np.linalg.solve(D, np.column_stack([b, A @ x0[::-1]]))
            elif op == "solveT":
                v = np.linalg.solve(D.T, b)
            elif op == "det":
                v = np.linalg.det(D)
            elif op == "inv":
                v = np.linalg.inv(D)
            elif op == "rank":
                v = int(np.linalg.matrix_rank(D)) if D.size else 0
            elif op == "nullcount":
                v = m - (int(np.linalg.matrix_rank(D)) if D.size else 0)
            elif op == "lstsq":
                v = np.linalg.lstsq(D, b, rcond=None)[0]
            elif op == "pinv":
                v = np.linalg.pinv(D)
            elif op in ("eig", "eigk", "eigmk"):
                if large:
                    which = "LM"
                    if j["sym"]:
                        ev = (spla.eigsh(A.tocsc(), k=k, which=which, return_eigenvectors=False)
                              if op == "eigk" else
                              spla.eigsh(A.tocsc(), k=k, sigma=0, which=which,
                                         return_eigenvectors=False))
                    else:
                        ev = (spla.eigs(A.tocsc(), k=k, which=which, return_eigenvectors=False)
                              if op == "eigk" else
                              spla.eigs(A.tocsc(), k=k, sigma=0, which=which,
                                        return_eigenvectors=False))
                    v = np.asarray(ev)
                else:
                    ev = np.linalg.eigvalsh(D) if j["sym"] else np.linalg.eigvals(D)
                    ev = ev[np.argsort(-np.abs(ev), kind="stable")]
                    v = ev if op == "eig" else (ev[:k] if op == "eigk" else ev[::-1][:k])
            elif op in ("svl", "svlk", "svdall", "svdk"):
                if large:
                    v = np.sort(spla.svds(A.tocsc(), k=k, return_singular_vectors=False))[::-1]
                else:
                    s = np.linalg.svd(D, compute_uv=False) if D.size else np.zeros(0)
                    if op in ("svdall", "svdk"):
                        v = s if op == "svdall" else s[:k]
                    else:
                        s = s[s > (s[0] if s.size else 0) * 1e-12] if s.size else s
                        v = s if op == "svl" else s[:k]
            elif op == "chol":
                v = sla.cholesky(D, lower=False)
            elif op == "charpoly":
                n = D.shape[0]
                v = (np.poly(D)[::-1] * (-1) ** n)
            elif op == "expm":
                v = sla.expm(D)
            elif op == "expmv":
                v = spla.expm_multiply(A.tocsc(), b)
            elif op == "logm":
                v = sla.logm(D)
            elif op == "sinm":
                v = sla.sinm(D)
            elif op == "pow2":
                v = D @ D
            elif op == "sqrtm":
                v = sla.sqrtm(D)
            elif op == "adj":
                v = np.linalg.det(D) * np.linalg.inv(D)
            elif op in ("norm1", "norm2", "normInf", "normF"):
                o = {"norm1": 1, "norm2": 2, "normInf": np.inf, "normF": "fro"}[op]
                if large:
                    v = spla.norm(A, o)
                else:
                    v = np.linalg.norm(D, o) if D.size else 0.0
            elif op == "tr":
                v = A.diagonal().sum()
            elif op == "dotAAt":
                v = D @ D.T
            elif op == "dotv":
                v = A @ x0
            elif op == "T":
                v = D.T
            elif op == "kron":
                v = np.kron(D, K)
            elif op == "symq":
                v = bool(np.array_equal(D, D.T))
            elif op == "hermq":
                v = bool(np.array_equal(D, D.conj().T))
            elif op == "pdq":
                v = bool(np.all(np.linalg.eigvalsh((D + D.conj().T) / 2) > 0))
            if v is not None:
                if (isinstance(v, (np.ndarray, float, complex, np.floating, np.complexfloating))
                        and not np.all(np.isfinite(v))):
                    out["err"] = "non-finite"
                else:
                    out["val"] = enc(v)
        except Exception as e:   # noqa: BLE001 -- one bad case must not stop the batch
            out["err"] = f"{type(e).__name__}: {e}"[:200]
        print("@@S " + json.dumps(out), flush=True)


# --------------------------------------------------------------------------
# maxima
# --------------------------------------------------------------------------

def maxima_cases(cases):
    return [c for c in cases if c.op.maxima and c.group in ("exact", "sym") and c.op.unique(c)
            and c.rows <= 8 and c.cols <= 8
            and not (c.op.maxima == "eig" and (c.rows > 3 or c.kind == "sym"))]


def run_maxima(cases, workdir, maxima, log):
    sel = maxima_cases(cases)
    if not sel:
        return {}
    path = os.path.join(workdir, "maxima_cases.mac")
    with open(path, "w") as fh:
        for c in sel:
            fh.write(f'SLrun({c.id}, "{c.op.maxima}", {maxima_matrix(c.mat)}, '
                     f'{"true" if c.kind == "sym" else "false"})$\n')
    t0 = time.time()
    batch = f'batchload("{MAXIMA_DRIVER}")$ batchload("{path}")$'
    rc, out, err = _run([maxima, "--very-quiet", f"--batch-string={batch}"], 120 + 20 * len(sel))
    with open(os.path.join(workdir, "maxima.out"), "w") as fh:
        fh.write(out + err)
    res = parse_result_lines(out, "@@M ")
    log(f"maxima: {len(res)}/{len(sel)} answers in {time.time() - t0:.0f} s")
    return res


# --------------------------------------------------------------------------
# Values: parsing and comparison
# --------------------------------------------------------------------------

getcontext().prec = 80
NUM_RE = re.compile(r"^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)(`(`?)([0-9.]*))?"
                    r"(?:\*\^([-+]?\d+))?\s*$")


class Num:
    """A number: re, im are Fraction (exact / high precision) or float."""
    __slots__ = ("re", "im", "kind", "prec")

    def __init__(self, re_, im, kind, prec=None):
        self.re, self.im, self.kind, self.prec = re_, im, kind, prec

    def __repr__(self):
        def f(x):
            return f"{float(x):.10g}" if not isinstance(x, Fraction) or x.denominator != 1 else str(x)
        return f(self.re) if self.im == 0 else f"{f(self.re)}{'+' if self.im >= 0 else ''}{f(self.im)}I"


class Sym:
    __slots__ = ("s",)

    def __init__(self, s):
        self.s = s

    def __repr__(self):
        return self.s[:60]


def pnum(s):
    if not isinstance(s, str) or len(s) < 2 or s[1] != ":":
        return Sym(str(s))
    tag, body = s[0], s[2:]
    try:
        if tag == "i":
            return Num(Fraction(int(body)), Fraction(0), "e")
        if tag == "q":
            return Num(Fraction(body), Fraction(0), "e")
        if tag in ("r", "R"):
            if body in ("nan", "inf", "-inf"):
                return Sym(body)
            mt = NUM_RE.match(body)
            if not mt:
                return Sym(body)
            mant, tick, acc, pr, ex = mt.groups()
            txt = mant + ("e" + ex if ex else "")
            if tag == "r":
                return Num(float(txt), 0.0, "m")
            prec = float(pr) if (tick and not acc and pr) else None
            return Num(Fraction(txt), Fraction(0), "h", prec)
    except (ValueError, ZeroDivisionError):
        pass
    return Sym(body)


def parse_val(v):
    if v is None:
        return None
    if isinstance(v, bool):
        return v
    if isinstance(v, list):
        return [parse_val(x) for x in v]
    if isinstance(v, dict) and "c" in v:
        re_, im = parse_val(v["c"][0]), parse_val(v["c"][1])
        if not isinstance(re_, Num) or not isinstance(im, Num):
            return Sym(json.dumps(v))
        kind = "m" if "m" in (re_.kind, im.kind) else ("h" if "h" in (re_.kind, im.kind) else "e")
        precs = [p for p in (re_.prec, im.prec) if p]
        return Num(re_.re, im.re, kind, min(precs) if precs else None)
    if isinstance(v, (int, float)):
        return Num(float(v), 0.0, "m")
    return pnum(v)


def leaves(t):
    if isinstance(t, list):
        for x in t:
            yield from leaves(x)
    else:
        yield t


def d2(a, b):
    """|a - b|^2 as a float (exact differences for exact / high precision)."""
    if a.kind == "m" or b.kind == "m":
        dr = float(a.re) - float(b.re)
        di = float(a.im) - float(b.im)
        return dr * dr + di * di
    try:
        return float((a.re - b.re) ** 2 + (a.im - b.im) ** 2)
    except OverflowError:
        return float("inf")


def m2(a):
    try:
        if a.kind == "m":
            return float(a.re) ** 2 + float(a.im) ** 2
        return float(a.re * a.re + a.im * a.im)
    except OverflowError:
        return float("inf")


def nabs(a):
    """|a| as a Num (Decimal sqrt for exact / high precision)."""
    if a.im == 0:
        return Num(abs(a.re), a.im * 0, a.kind, a.prec)
    if a.kind == "m":
        return Num(math.hypot(a.re, a.im), 0.0, "m")
    x = a.re * a.re + a.im * a.im
    r = (Decimal(x.numerator) / Decimal(x.denominator)).sqrt()
    return Num(Fraction(r), Fraction(0), "h", a.prec)


class Tol:
    def __init__(self, rtol, atol=0.0, exact=False, mpfr=False):
        self.rtol, self.atol, self.exact, self.mpfr = rtol, atol, exact, mpfr


def cmp_tree(a, b, tol):
    """Normwise: max|a_i - b_i| <= rtol max|b_i| + atol, exact leaves equal
    when tol.exact.  Returns (ok, detail)."""
    pairs = []

    def walk(x, y, path):
        if isinstance(y, list):
            if not isinstance(x, list):
                raise _Mismatch(f"shape: got {short(x)} at {path}, expected a list of {len(y)}")
            if len(x) != len(y):
                raise _Mismatch(f"shape: length {len(x)} vs {len(y)} at {path or 'top'}")
            for i, (p, q) in enumerate(zip(x, y)):
                walk(p, q, path + f"[{i + 1}]")
        elif isinstance(y, bool) or isinstance(x, bool):
            if x is not y:
                raise _Mismatch(f"{short(x)} vs {short(y)} at {path or 'top'}")
        elif isinstance(y, Sym) or isinstance(x, Sym):
            if not (isinstance(x, Sym) and isinstance(y, Sym) and x.s == y.s):
                raise _Mismatch(f"{short(x)} vs {short(y)} at {path or 'top'}")
        elif isinstance(x, list):
            raise _Mismatch(f"shape: got a list at {path or 'top'}, expected a scalar")
        else:
            pairs.append((x, y, path))

    try:
        walk(a, b, "")
    except _Mismatch as e:
        return False, str(e)
    if tol.exact:
        for x, y, path in pairs:
            if x.kind == "m" and y.kind != "m":
                return False, f"machine number {x} in an exact result at {path or 'top'}"
            if x.kind == "e" and y.kind == "e" and (x.re != y.re or x.im != y.im):
                return False, f"{x} vs {y} at {path or 'top'} (exact)"
    if tol.mpfr:
        for x, y, path in pairs:
            if x.kind == "m" and y.kind == "h":
                return False, (f"machine number {x} in an MPFR result at {path or 'top'} "
                               f"(precision dropped to machine)")
    worst, where, scale2 = 0.0, None, 0.0
    for x, y, path in pairs:
        dd = d2(x, y)
        scale2 = max(scale2, m2(y))
        if not dd <= worst:
            worst, where = dd, (x, y, path)
    err, scale = math.sqrt(worst), math.sqrt(scale2)
    lim = tol.rtol * scale + tol.atol
    if err <= lim or (err == 0):
        return True, ""
    x, y, path = where
    rel = err / scale if scale else err
    return False, f"max err {err:.3g} (rel {rel:.3g}) > tol {lim:.3g}: {x} vs {y} at {path or 'top'}"


class _Mismatch(Exception):
    pass


def cmp_multiset(a, b, tol, absolute=False):
    if not isinstance(a, list) or not isinstance(b, list):
        return cmp_tree(a, b, tol)
    xa = [x for x in leaves(a)]
    xb = [x for x in leaves(b)]
    if any(not isinstance(x, Num) for x in xa + xb):
        return cmp_tree(a, b, tol)
    if absolute:
        xa, xb = [nabs(x) for x in xa], [nabs(x) for x in xb]
    if len(xa) != len(xb):
        return False, f"{len(xa)} values vs {len(xb)}: {short(xa)} vs {short(xb)}"
    if tol.exact and all(x.kind == "e" for x in xa + xb):
        key = lambda x: (x.re, x.im)   # noqa: E731
        ok = sorted(xa, key=key) == sorted(xb, key=key) if False else \
            sorted((x.re, x.im) for x in xa) == sorted((x.re, x.im) for x in xb)
        return (True, "") if ok else (False, f"{short(xa)} vs {short(xb)} (exact multiset)")
    scale = math.sqrt(max([m2(y) for y in xb] + [0.0]))
    lim = tol.rtol * scale + tol.atol
    left = list(range(len(xb)))
    for x in sorted(xa, key=lambda z: -m2(z)):
        best, bi = None, None
        for idx in left:
            dd = d2(x, xb[idx])
            if best is None or dd < best:
                best, bi = dd, idx
        if best is None or math.sqrt(best) > lim:
            return False, (f"{x} has no partner within {lim:.3g} (nearest "
                           f"{xb[bi] if bi is not None else '-'}): {short(xa)} vs {short(xb)}")
        left.remove(bi)
    return True, ""


def short(x, n=160):
    s = repr(x)
    return s if len(s) <= n else s[:n] + "..."


def min_prec(t):
    ps = [x.prec for x in leaves(t) if isinstance(x, Num) and x.prec]
    return min(ps) if ps else None


def has_nonfinite(t):
    return any(isinstance(x, Sym) and ("inf" in x.s or "nan" in x.s) for x in leaves(t))


# --------------------------------------------------------------------------
# Verdicts
# --------------------------------------------------------------------------

ORDER = ["CRASH", "TIMEOUT", "FAIL", "UNEVAL", "PASS", "SKIP"]
RANK = {v: i for i, v in enumerate(ORDER)}


def cond_of(case):
    if case.cond is not None and not (isinstance(case.cond, float) and math.isnan(case.cond)):
        return case.cond
    return DEFAULT_COND.get(case.cls, 1e3)


def value_tol(case, oracle_vals):
    """Tolerance for a value comparison against one oracle's answer."""
    op = case.op
    cond = cond_of(case) if op.cs else 1.0
    cond = max(1.0, cond if math.isfinite(cond) else 1e12)
    oracle_machine = any(isinstance(x, Num) and x.kind == "m" for x in leaves(oracle_vals))
    if op.id in ITERATIVE:
        # a Krylov answer is as accurate as its stopping tolerance, not the
        # working precision: MPFR input does not make it 25-digit
        oracle_machine = True
    if case.group in ("exact", "sym") and not oracle_machine:
        return Tol(1e-30, exact=True)
    if case.group == "mpfr" and not oracle_machine:
        p = min_prec(oracle_vals)
        base = min(1e-25 * op.tol * cond, 1e-10)
        if p:
            base = max(base, 10.0 ** -(p - 2))
        return Tol(base, mpfr=True)
    return Tol(min(1e-8 * op.tol * cond, 1e-3 * op.tol))


def scipy_tol(case):
    op = case.op
    cond = cond_of(case) if op.cs else 1.0
    cond = max(1.0, cond if math.isfinite(cond) else 1e12)
    return Tol(min(1e-8 * op.tol * cond, 1e-3 * op.tol))


def inv_tol(case):
    op = case.op
    cond = max(1.0, cond_of(case)) if op.inv_cond else 1.0
    if not math.isfinite(cond):
        cond = 1e12
    if case.group in ("exact", "sym"):
        return 1e-30
    if case.group == "mpfr":
        return min(1e-20 * op.tol * cond, 1e-8)
    return min(1e-8 * op.tol * cond, 1e-3)


def compare(case, mine, theirs, tol, mode=None):
    mode = mode or case.op.cmp
    if mode == "set":
        return cmp_multiset(mine, theirs, tol)
    if mode == "absset":
        return cmp_multiset(mine, theirs, tol, absolute=True)
    if mode in ("count", "bool"):
        return cmp_tree(mine, theirs, Tol(0.0, exact=True))
    if mode == "det" and not tol.exact:
        return cmp_tree(mine, theirs, Tol(tol.rtol, 1e-12 * case.op.tol * hadamard(case.mat)))
    return cmp_tree(mine, theirs, tol)


def msg_tags(msgs):
    return [m for m in (msgs or []) if m not in ("General::stop",)]


def evaluate_case(case, mt, crash, ref, sp, mx, leak):
    """Rows (oracle, verdict, detail) for one case."""
    op = case.op
    rows = []

    def row(oracle, verdict, detail=""):
        if verdict == "FAIL":
            why = exempt_reason(op.id, case.cls, case.kind, oracle)
            if why:
                verdict, detail = "SKIP", f"EXEMPT: {why} | {detail}"
        rows.append((oracle, verdict, detail))

    if leak:
        row("valgrind", "FAIL", leak)
    if crash:
        row("mathilda", crash[0], crash[1])
        return rows
    if mt is None:
        row("mathilda", "CRASH", "no answer from Mathilda (lost output)")
        return rows
    sig = mt["sig"]
    if sig == "T":
        row("mathilda", "TIMEOUT", f"> per-case TimeConstrained")
        return rows
    rsig = ref.get("sig") if ref else None
    oracle_answers = (rsig not in (None, "U", "E", "T")) or (sp and "val" in sp) or (mx and "val" in mx)
    my_msgs, ref_msgs = msg_tags(mt.get("msgs")), msg_tags(ref.get("msgs") if ref else None)

    def msg_row():
        if ref is None:
            return
        want = sorted({m.split("::")[-1] for m in ref_msgs})
        have = {m.split("::")[-1] for m in my_msgs}
        missing = [t for t in want if t not in have]
        if missing:
            row("msg", "FAIL", f"reference emits {ref_msgs}; Mathilda emits {my_msgs or 'none'}")
        elif want or my_msgs:
            row("msg", "PASS", f"ref {ref_msgs} mathilda {my_msgs}")
        else:
            row("msg", "PASS")

    if sig == "U":
        if rsig in ("U", "E"):
            # the reference defines the semantics: it does not evaluate this
            # call either (scipy answering a different question is no evidence)
            row("sig", "PASS" if rsig == "U" else "FAIL",
                "both unevaluated" if rsig == "U" else "reference gives $Failed, Mathilda U")
            msg_row()
        elif oracle_answers:
            row("mathilda", "UNEVAL", f"reference sig {rsig}" if rsig else "an oracle answers it")
        else:
            row("mathilda", "UNEVAL", "no oracle to tell")
        return rows

    # Mathilda answered
    if ref is not None and rsig != "T":
        if rsig == sig:
            row("sig", "PASS", sig)
        else:
            row("sig", "FAIL", f"Mathilda {sig}, reference {rsig}")
        msg_row()
    mine = parse_val(mt.get("val"))
    if isinstance(mine, Sym) and "$TIMEOUT" in mine.s:
        row("mathilda", "TIMEOUT", "encoding the answer timed out")
        return rows
    uniq = op.cmp is not None and op.unique(case)

    # when the reference and scipy disagree with EACH OTHER the problem is
    # beyond machine precision (or the two define it differently): neither is
    # an oracle for this case
    ref_v = parse_val(ref["val"]) if (uniq and ref is not None and rsig not in ("U", "E", "T")
                                      and ref.get("val") is not None) else None
    sp_v = parse_val(sp["val"]) if (uniq and sp is not None and sp.get("val") is not None) else None
    if sp_v is not None and has_nonfinite(sp_v):
        row("scipy", "SKIP", "scipy answer not finite")
        sp_v = None
    disagree = ""
    if ref_v is not None and sp_v is not None and case.group == "machine":
        ok, det = compare(case, ref_v, sp_v, scipy_tol(case))
        if not ok:
            disagree = "reference and scipy disagree with each other (" + det + ")"
    if ref_v is not None:
        ok, det = compare(case, mine, ref_v, value_tol(case, ref_v))
        if disagree and not ok:
            row("ref", "SKIP", disagree)
        else:
            row("ref", "PASS" if ok else "FAIL", det)
    if sp_v is not None:
        ok, det = compare(case, mine, sp_v, scipy_tol(case))
        if disagree and not ok:
            row("scipy", "SKIP", disagree)
        else:
            row("scipy", "PASS" if ok else "FAIL", det)
    if uniq and mx is not None and mx.get("val") is not None:
        theirs = parse_val(mx["val"])
        if any(isinstance(x, Sym) for x in leaves(theirs)):
            row("maxima", "SKIP", "maxima answer not numeric: " + short(theirs))
        else:
            ok, det = compare(case, mine, theirs, value_tol(case, theirs))
            row("maxima", "PASS" if ok else "FAIL", det)

    if op.inv:
        iv = parse_val(mt.get("inv"))
        tol = inv_tol(case)
        ref_iv = parse_val(ref.get("inv")) if ref else None
        ref_ok = isinstance(ref_iv, Num) and math.sqrt(m2(ref_iv)) <= tol
        if isinstance(iv, Sym) and "$TIMEOUT" in iv.s:
            row("inv", "TIMEOUT", "computing the invariant (Mathilda's own arithmetic on the "
                                  "answer) exceeded the per-case limit")
        elif isinstance(iv, Num) and math.sqrt(m2(iv)) <= tol:
            row("inv", "PASS", f"residual {float(iv.re):.3g}")
        elif ref is not None and rsig not in ("U", "E", "T") and not ref_ok:
            row("inv", "SKIP", f"invariant does not apply: the reference fails it too "
                               f"(mathilda {short(iv, 80)}, reference {short(ref_iv, 80)})")
        elif isinstance(iv, Sym) and "MatrixExp[" in iv.s:
            row("inv", "SKIP", "invariant needs MatrixExp, unevaluated")
        else:
            row("inv", "FAIL", f"residual {short(iv, 200)} > {tol:.3g}")

    if op.order and isinstance(mine, list):
        xs = [x for x in leaves(mine) if isinstance(x, Num)]
        mags = [math.sqrt(m2(x)) for x in xs]
        lim = 1e-8 * max(mags + [0.0]) * (1 if case.group != "mpfr" else 1e-10)
        if case.group in ("exact", "sym"):
            lim = 1e-30 * max(mags + [0.0])
        bad = [i for i in range(len(mags) - 1) if mags[i + 1] > mags[i] + lim]
        row("order", "FAIL" if bad else "PASS",
            f"|values| increase at position {bad[0] + 1}: {short(xs)}" if bad else "")

    if op.self_ and op.cmp and case.tier != "large" and mt.get("dsig") is not None:
        dsig = mt.get("dsig")
        if dsig in ("U", "E", "T") or mt.get("dval") is None:
            row("self", "SKIP", f"dense call: {dsig}")
        else:
            dense = parse_val(mt["dval"])
            ok, det = compare(case, mine, dense, value_tol(case, dense))
            if ok:
                row("self", "PASS")
            else:
                ref_differs = False
                if ref and ref.get("val") is not None and ref.get("dval") is not None:
                    rv, rd = parse_val(ref["val"]), parse_val(ref["dval"])
                    ref_differs = not compare(case, rv, rd, value_tol(case, rd))[0]
                if ref_differs:
                    row("self", "PASS", "reference also answers sparse and dense differently: " + det)
                else:
                    row("self", "FAIL", "Normal[f[S]] != f[Normal[S]]: " + det)
    if not any(r[0] in ("ref", "scipy", "maxima", "inv", "self", "sig") for r in rows):
        row("value", "SKIP", "no oracle answer to compare")
    return rows


# --------------------------------------------------------------------------
# Recording (--record / --quick)
# --------------------------------------------------------------------------

def save_record(cases, ref, sp, mx):
    os.makedirs(RECORD_DIR, exist_ok=True)
    data = {"meta": {"seed": QUICK_SEED, "n": QUICK_N, "gen_version": GEN_VERSION,
                     "recorded": time.strftime("%Y-%m-%d"),
                     "note": "oracle answers for tools/sparse_linalg_oracle.py --quick; "
                             "regenerate with --record"},
            "cases": {}}
    for c in cases:
        e = {"op": c.op.id, "cls": c.cls, "kind": c.kind, "size": c.size}
        if c.id in ref:
            e["ref"] = {k: ref[c.id].get(k) for k in ("sig", "msgs", "val", "inv", "dsig", "dval")}
        if c.id in sp:
            e["scipy"] = {k: sp[c.id][k] for k in ("val", "cond", "err") if k in sp[c.id]}
        if c.id in mx:
            e["maxima"] = {k: mx[c.id][k] for k in ("val", "err") if k in mx[c.id]}
        data["cases"][c.fp] = e
    with open(RECORD_FILE, "w") as fh:
        json.dump(data, fh, separators=(",", ":"), sort_keys=True)
        fh.write("\n")
    return RECORD_FILE


def load_record(cases):
    with open(RECORD_FILE) as fh:
        data = json.load(fh)
    ref, sp, mx = {}, {}, {}
    missing = 0
    for c in cases:
        e = data["cases"].get(c.fp)
        if e is None:
            missing += 1
            continue
        if "ref" in e:
            ref[c.id] = dict(e["ref"], id=c.id)
        if "scipy" in e:
            sp[c.id] = dict(e["scipy"], id=c.id)
        if "maxima" in e:
            mx[c.id] = dict(e["maxima"], id=c.id)
    return ref, sp, mx, missing


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def parse_list(s):
    return [x.strip() for x in s.split(",") if x.strip()] if s else None


def repro_text(c):
    return (f"(* case {c.id}: op {c.op.id}, class {c.cls}, kind {c.kind}, {c.size} *)\n"
            f"A = {c.wl};\n"
            f"x0 = Table[(-1)^i i/(i + 1), {{i, {c.cols}}}]"
            + {"machine": " // N", "mpfr": " // N[#, 30] &"}.get(c.group, "") + ";\n"
            f"b = Normal[A . x0]; Bm = Transpose[{{b, Normal[A . Reverse[x0]]}}];\n"
            f"{c.tpl}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", type=int, default=300, help="number of small/medium cases (default 300)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--large", type=int, nargs="?", const=6, default=0,
                    help="add N large cases (2000..20000, machine kinds, scipy oracle; default 6)")
    ap.add_argument("--medium-frac", type=float, default=0.2)
    ap.add_argument("--only", help="comma-separated op ids (globs ok)")
    ap.add_argument("--classes", help="comma-separated classes: " + ",".join(CLASSES))
    ap.add_argument("--kinds", help="comma-separated kinds: " + ",".join(KINDS))
    ap.add_argument("--oracles", default="ref,scipy,maxima")
    ap.add_argument("--quick", action="store_true", help="replay recorded oracle answers")
    ap.add_argument("--record", action="store_true", help="record oracle answers for --quick")
    ap.add_argument("--tsv", help="TSV report path (default: <workdir>/report.tsv)")
    ap.add_argument("--keep", metavar="DIR", help="keep the generated programs and MTX files in DIR")
    ap.add_argument("--no-fail", action="store_true")
    ap.add_argument("--chunk", type=int, default=100)
    ap.add_argument("-j", "--jobs", type=int, default=4, help="parallel Mathilda processes")
    ap.add_argument("--timeout", type=float, default=10, help="Mathilda per-case TimeConstrained")
    ap.add_argument("--ref-timeout", type=float, default=30)
    ap.add_argument("--binary", default=BINARY)
    ap.add_argument("--wolframscript", default=WOLFRAMSCRIPT)
    ap.add_argument("--python", default=PYTHON)
    ap.add_argument("--maxima", default=MAXIMA)
    ap.add_argument("--valgrind", action="store_true")
    ap.add_argument("--repro", type=int, metavar="ID", help="print the WL input of case ID and exit")
    ap.add_argument("--show", type=int, default=3,
                    help="print up to N repros per failing (op, oracle) group (default 3)")
    ap.add_argument("--list-ops", action="store_true")
    ap.add_argument("--heads", action="store_true",
                    help="check Mathilda's result heads against the checked-in reference table "
                         "tools/data/sparse_linalg_heads.tsv (runs sparse_linalg_semantics.py)")
    ap.add_argument("--scipy-worker", help=argparse.SUPPRESS)
    a = ap.parse_args()

    if a.scipy_worker:
        scipy_worker(a.scipy_worker)
        return 0
    if a.heads:
        table = os.path.join(ROOT, "tools", "data", "sparse_linalg_heads.tsv")
        if not os.path.exists(table):
            sys.exit(f"no {table}: generate it with tools/sparse_linalg_semantics.py")
        return subprocess.call([sys.executable, os.path.join(ROOT, "tools", "sparse_linalg_semantics.py"),
                                "--mathilda", a.binary, "--no-reference", "--out", table])
    if a.list_ops:
        for o in OPS:
            print(f"{o.id:32s} {o.tpl}")
        return 0

    ops = OPS
    if a.only:
        pats = parse_list(a.only)
        ops = [o for o in OPS if any(fnmatch.fnmatch(o.id, p) for p in pats)]
        if not ops:
            sys.exit("no op matches --only")
    if a.quick or a.record:
        a.seed, a.n, a.large, a.medium_frac = QUICK_SEED, QUICK_N, 0, 0.0
        cases = gen_cases(QUICK_N, QUICK_SEED, medium_frac=0.0)
    else:
        cases = gen_cases(a.n, a.seed, a.medium_frac, a.large, ops,
                          parse_list(a.classes), parse_list(a.kinds))
    if a.repro is not None:
        for c in cases:
            if c.id == a.repro:
                print(repro_text(c))
                return 0
        sys.exit(f"no case {a.repro}")

    if not os.path.exists(a.binary):
        # exit 77: ctest's SKIP_RETURN_CODE (the CMake test tree builds no REPL;
        # the oracle drives the top-level ./Mathilda from `make`)
        print(f"SKIP: no Mathilda binary at {a.binary} (run `make` at the top level)")
        sys.exit(77)
    workdir = os.path.abspath(a.keep) if a.keep else tempfile.mkdtemp(prefix="sla_oracle_")
    os.makedirs(workdir, exist_ok=True)
    t_start = time.time()

    def log(msg):
        sys.stderr.write(f"[{time.time() - t_start:6.0f}s] {msg}\n")
        sys.stderr.flush()

    log(f"{len(cases)} cases (seed {a.seed}); workdir {workdir}")
    oracles = set(parse_list(a.oracles) or [])
    ref, sp, mx = {}, {}, {}
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = {}
        if a.quick:
            if not os.path.exists(RECORD_FILE):
                sys.exit(f"no recording at {RECORD_FILE}; run --record")
            ref, sp, mx, missing = load_record(cases)
            if missing:
                log(f"WARNING: {missing} of {len(cases)} cases are not in the recording "
                    f"(generator changed?) -- re-run --record")
                if missing > len(cases) // 10:
                    print(f"stale recording: {missing}/{len(cases)} cases missing; run --record")
                    return 2
        else:
            small = [c for c in cases if c.tier != "large"]
            if "ref" in oracles and os.path.exists(a.wolframscript):
                futs["ref"] = ex.submit(run_reference, small, workdir, a.ref_timeout,
                                        a.wolframscript, log)
            elif "ref" in oracles:
                log(f"reference: no wolframscript at {a.wolframscript}; skipped")
            if "scipy" in oracles and os.path.exists(a.python):
                futs["scipy"] = ex.submit(run_scipy, cases, workdir, a.python, log)
            elif "scipy" in oracles:
                log(f"scipy: no interpreter at {a.python}; skipped")
            if "maxima" in oracles and a.maxima and os.path.exists(a.maxima):
                futs["maxima"] = ex.submit(run_maxima, cases, workdir, a.maxima, log)
        runner = MathildaRunner(a.binary, workdir, a.timeout, a.valgrind, a.keep)
        chunk = max(1, a.chunk // 10) if a.valgrind else a.chunk
        mres, crashes, leaks = runner.run_all(cases, chunk, a.jobs)
        log(f"mathilda: {len(mres)}/{len(cases)} answers, {len(crashes)} crash/timeout")
        if "ref" in futs:
            ref = futs["ref"].result()
        if "scipy" in futs:
            sp = futs["scipy"].result()
        if "maxima" in futs:
            mx = futs["maxima"].result()

    for c in cases:
        s = sp.get(c.id)
        if s and "cond" in s:
            try:
                c.cond = float(s["cond"])
            except (TypeError, ValueError):
                pass

    if a.record:
        path = save_record(cases, ref, sp, mx)
        log(f"recorded {len(ref)} reference, {len(sp)} scipy, {len(mx)} maxima answers -> {path}")

    # ---- verdicts -------------------------------------------------------
    tsv_path = a.tsv or os.path.join(workdir, "report.tsv")
    all_rows = []
    case_verdict = {}
    for c in cases:
        rows = evaluate_case(c, mres.get(c.id), crashes.get(c.id), ref.get(c.id),
                             sp.get(c.id), mx.get(c.id), leaks.get(c.id))
        worst = min((r[1] for r in rows), key=lambda v: RANK[v], default="SKIP")
        case_verdict[c.id] = worst
        for oracle, verdict, detail in rows:
            all_rows.append((c, oracle, verdict, detail))
    with open(tsv_path, "w") as fh:
        fh.write("case\tclass\tkind\tsize\top\toracle\tverdict\tdetail\n")
        for c, oracle, verdict, detail in all_rows:
            fh.write(f"{c.id}\t{c.cls}\t{c.kind}\t{c.size}\t{c.op.id}\t{oracle}\t{verdict}\t"
                     f"{detail.replace(chr(9), ' ').replace(chr(10), ' ')[:600]}\n")

    # ---- summary ---------------------------------------------------------
    from collections import Counter, defaultdict
    per_op = defaultdict(Counter)
    for c in cases:
        per_op[c.op.id][case_verdict[c.id]] += 1
    cols = ORDER
    print(f"\nper op (case verdict = worst row), {len(cases)} cases, seed {a.seed}:")
    print(f"{'op':32s}" + "".join(f"{v:>8s}" for v in cols))
    tot = Counter()
    for oid in [o.id for o in OPS if o.id in per_op]:
        cnt = per_op[oid]
        tot.update(cnt)
        print(f"{oid:32s}" + "".join(f"{cnt.get(v, 0) or '.':>8}" for v in cols))
    print(f"{'TOTAL':32s}" + "".join(f"{tot.get(v, 0):>8}" for v in cols))
    per_or = defaultdict(Counter)
    for c, oracle, verdict, _ in all_rows:
        per_or[oracle][verdict] += 1
    print("\nper check (rows):")
    print(f"{'check':32s}" + "".join(f"{v:>8s}" for v in cols))
    for oracle in sorted(per_or):
        print(f"{oracle:32s}" + "".join(f"{per_or[oracle].get(v, 0) or '.':>8}" for v in cols))

    groups = defaultdict(list)
    for c, oracle, verdict, detail in all_rows:
        if verdict in ("FAIL", "CRASH", "TIMEOUT"):
            groups[(c.op.id, oracle, verdict)].append((c, detail))
    if groups and a.show:
        print("\nfailures (first per op x check; --repro ID prints a case):")
        for (oid, oracle, verdict), items in sorted(groups.items()):
            print(f"\n== {verdict} {oid} [{oracle}] x{len(items)}")
            for c, detail in items[:a.show]:
                wl = c.wl if len(c.wl) < 300 else c.wl[:300] + "..."
                print(f"  #{c.id} {c.cls}/{c.kind} {c.size}: {detail[:400]}")
                print(f"      {c.tpl}   with A = {wl}")
    print(f"\nreport: {tsv_path}   ({time.time() - t_start:.0f} s)")
    bad = tot.get("FAIL", 0) + tot.get("CRASH", 0)
    return 1 if bad and not a.no_fail else 0


if __name__ == "__main__":
    sys.exit(main())
