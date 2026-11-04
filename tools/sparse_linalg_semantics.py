#!/usr/bin/env python3
"""
sparse_linalg_semantics.py -- what does each linear-algebra head return when
its argument is a SparseArray?

The reference system is not consistent about this. Inverse returns a List;
MatrixExp and Normalize return a SparseArray; SingularValueDecomposition returns
{List, SparseArray, List}; Inverse of a symbolic SparseArray is left
unevaluated. The answer can also change with size or density. This tool asks
the reference instead of guessing. Every head template is run on matrices of
several entry kinds (integer, real, complex, 30-digit real, symbolic) and three
sizes (4, 40, 400), and the *shape signature* of each answer is recorded:

    S          a SparseArray
    L          a List that is a full array (vector, matrix, ...)
    x          a scalar (number, symbol, polynomial)
    U          the call came back unevaluated (same head as the call)
    E          an error / $Failed / $Aborted / timeout
    {a,b,..}   a ragged List whose parts have signatures a, b, ...
    Head       any other object (LinearSolveFunction, ...)

The table is written to tools/data/sparse_linalg_heads.tsv (one row per
id, kind and size). Mathilda's native sparse handlers wrap their answers to
match it, and tools/sparse_linalg_oracle.py --heads checks them against it.
Pass --mathilda to run the same probe through Mathilda and print the rows that
differ.

Usage:
    python3 tools/sparse_linalg_semantics.py [--wolframscript PATH]
        [--mathilda ./Mathilda] [--out tools/data/sparse_linalg_heads.tsv]
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (id, template). In a template, A is the square test matrix (nonsingular,
# diagonally dominant), P its symmetric positive-definite twin, T a 3n x n
# tall matrix, b a dense vector, B a sparse vector and n the size.
TEMPLATES = [
    ("Det", "Det[A]"),
    ("Inverse", "Inverse[A]"),
    ("LinearSolve-v", "LinearSolve[A, b]"),
    ("LinearSolve-sv", "LinearSolve[A, B]"),
    ("LinearSolve-m", "LinearSolve[A, Normal[A]]"),
    ("LinearSolve-sm", "LinearSolve[A, A]"),
    ("LinearSolve-f", "LinearSolve[A]"),
    ("LinearSolve-fapply", "LinearSolve[A][b]"),
    ("LinearSolve-fapplym", "LinearSolve[A][A]"),
    ("LinearSolve-krylov", "LinearSolve[P, b, Method -> \"Krylov\"]"),
    ("LinearSolve-chol", "LinearSolve[P, b, Method -> \"Cholesky\"]"),
    ("LinearSolve-mf", "LinearSolve[A, b, Method -> \"Multifrontal\"]"),
    ("Eigenvalues", "Eigenvalues[P]"),
    ("Eigenvalues-k", "Eigenvalues[P, 2]"),
    ("Eigenvalues-mk", "Eigenvalues[P, -2]"),
    ("Eigenvectors", "Eigenvectors[P]"),
    ("Eigenvectors-k", "Eigenvectors[P, 2]"),
    ("Eigensystem", "Eigensystem[P]"),
    ("Eigensystem-k", "Eigensystem[P, 2]"),
    ("Eigenvalues-gen", "Eigenvalues[{A, P}]"),
    ("MatrixRank", "MatrixRank[A]"),
    ("NullSpace", "NullSpace[A - A]"),
    ("NullSpace-1", "NullSpace[SparseArray[{{1, 1} -> 1, {1, 2} -> 1}, {2, n}]]"),
    ("RowReduce", "RowReduce[A]"),
    ("LUDecomposition", "LUDecomposition[A]"),
    ("QRDecomposition", "QRDecomposition[A]"),
    ("QRDecomposition-tall", "QRDecomposition[T]"),
    ("SVD", "SingularValueDecomposition[A]"),
    ("SVD-k", "SingularValueDecomposition[A, 2]"),
    ("SingularValueList", "SingularValueList[A]"),
    ("SingularValueList-k", "SingularValueList[A, 2]"),
    ("Cholesky", "CholeskyDecomposition[P]"),
    ("Schur", "SchurDecomposition[A]"),
    ("Hessenberg", "HessenbergDecomposition[A]"),
    ("Jordan", "JordanDecomposition[P]"),
    ("CharacteristicPolynomial", "CharacteristicPolynomial[A, t]"),
    ("LeastSquares", "LeastSquares[T, Join[b, b, b]]"),
    ("PseudoInverse", "PseudoInverse[T]"),
    ("MatrixExp", "MatrixExp[A/Norm[A, 1]]"),
    ("MatrixExp-v", "MatrixExp[A/Norm[A, 1], b]"),
    ("MatrixLog", "MatrixLog[P]"),
    ("MatrixFunction", "MatrixFunction[Sin, P]"),
    ("MatrixPower-2", "MatrixPower[A, 2]"),
    ("MatrixPower-m1", "MatrixPower[A, -1]"),
    ("MatrixPower-half", "MatrixPower[P, 1/2]"),
    ("MatrixPower-v", "MatrixPower[A, 2, b]"),
    ("Minors", "Minors[A[[1 ;; 3, 1 ;; 3]]]"),
    ("Adjugate", "Adjugate[A[[1 ;; 3, 1 ;; 3]]]"),
    ("Orthogonalize", "Orthogonalize[A]"),
    ("Projection", "Projection[B, b]"),
    ("Normalize", "Normalize[B]"),
    ("Norm-2", "Norm[A]"),
    ("Norm-1", "Norm[A, 1]"),
    ("Norm-F", "Norm[A, \"Frobenius\"]"),
    ("Tr", "Tr[A]"),
    ("Dot", "A . A"),
    ("Dot-v", "A . b"),
    ("Transpose", "Transpose[A]"),
    ("ConjugateTranspose", "ConjugateTranspose[A]"),
    ("KroneckerProduct", "KroneckerProduct[A[[1 ;; 2, 1 ;; 2]], A[[1 ;; 2, 1 ;; 2]]]"),
    ("Cross", "Cross[B[[1 ;; 3]], Reverse[B[[1 ;; 3]]]]"),
    ("HermiteDecomposition", "HermiteDecomposition[Round[A]]"),
    ("PositiveDefiniteMatrixQ", "PositiveDefiniteMatrixQ[P]"),
    ("SymmetricMatrixQ", "SymmetricMatrixQ[P]"),
    ("HermitianMatrixQ", "HermitianMatrixQ[P]"),
    ("DiagonalMatrixQ", "DiagonalMatrixQ[A]"),
    ("UpperTriangularMatrixQ", "UpperTriangularMatrixQ[A]"),
    ("SquareMatrixQ", "SquareMatrixQ[A]"),
    ("TensorProduct", "TensorProduct[B, B]"),
    ("Diagonal", "Diagonal[A]"),
    ("DiagonalMatrix", "DiagonalMatrix[B]"),
    ("UpperTriangularize", "UpperTriangularize[A]"),
    ("Inner", "Inner[Times, A, b, Plus]"),
    ("Outer", "Outer[Times, B, B]"),
]

KINDS = ["int", "real", "complex", "mpfr", "sym"]
SIZES = [4, 40, 400]

PRELUDE = r"""
mk[kind_, n_] := Module[{rules, f},
  f = Switch[kind,
    "int", Function[{i, j}, Mod[7 i + 3 j, 5] - 2],
    "real", Function[{i, j}, N[Mod[7 i + 3 j, 5] - 2]/3],
    "complex", Function[{i, j}, N[(Mod[7 i + 3 j, 5] - 2) + I Mod[i + j, 3]]/3],
    "mpfr", Function[{i, j}, N[(Mod[7 i + 3 j, 5] - 2)/3, 30]],
    "sym", Function[{i, j}, If[i == 1 && j == 2, aa, Mod[7 i + 3 j, 5] - 2]]];
  rules = Join[
    Table[{i, i} -> 4 + f[i, i]^2 + If[kind === "sym", 0, 0], {i, n}],
    Table[{i, Mod[i, n] + 1} -> f[i, Mod[i, n] + 1], {i, n}],
    Table[{i, Mod[i + 2, n] + 1} -> f[i, i + 2], {i, n}]];
  SparseArray[rules, {n, n}]];
sig[r_] := Which[
  Head[r] === SparseArray, "S",
  r === $Failed || r === $Aborted || r === $TimedOut, "E",
  ListQ[r] && ArrayQ[r] && FreeQ[r, _SparseArray], "L",
  ListQ[r], "{" <> StringRiffle[sig /@ r, ","] <> "}",
  AtomQ[r] || NumericQ[r] || Head[r] === Plus || Head[r] === Times ||
    Head[r] === Power || Head[r] === Rational || Head[r] === Complex, "x",
  True, ToString[Head[r]]];
SetAttributes[R, HoldAll];
R[id_, kind_, nn_, hd_, e_] := Block[{A, P, T, b, B, n = nn, r},
  A = mk[kind, n]; P = (A + ConjugateTranspose[A])/2 + n IdentityMatrix[n, SparseArray];
  T = Join[A, A, A]; b = N[Range[n]]/n; If[kind === "int" || kind === "sym", b = Range[n]];
  If[kind === "mpfr", b = N[Range[n]/n, 30]];
  B = SparseArray[{1 -> 3, 3 -> 4}, n];
  r = TimeConstrained[Quiet[e], 30, $TimedOut];
  WriteString["stdout", id, "\t", kind, "\t", n, "\t",
    If[Head[r] === Symbol[hd] || Head[Head[r]] === Symbol[hd], "U", sig[r]], "\n"]];
"""


def program(templates, sizes, kinds):
    out = [PRELUDE]
    for tid, tpl in templates:
        for kind in kinds:
            for n in sizes:
                if kind in ("sym", "mpfr") and n > 40:
                    continue
                hd = "Dot" if " . " in tpl else tpl.split("[", 1)[0]
                out.append(f'R["{tid}", "{kind}", {n}, "{hd}", {tpl}];')
    return "\n".join(out) + "\n"


def run(cmd, path, timeout):
    try:
        p = subprocess.run(cmd + [path], capture_output=True, text=True, timeout=timeout)
        return p.stdout
    except subprocess.TimeoutExpired as e:
        return e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")


def parse(text):
    rows = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) == 4:
            rows[tuple(parts[:3])] = parts[3]
    return rows


def find_wolframscript(arg):
    for c in (arg, os.environ.get("WOLFRAMSCRIPT"), shutil.which("wolframscript"),
              "/Applications/Mathematica.app/Contents/MacOS/wolframscript"):
        if c and os.path.exists(c):
            return c
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wolframscript")
    ap.add_argument("--mathilda")
    ap.add_argument("--only", help="comma-separated template ids")
    ap.add_argument("--out", default=os.path.join(ROOT, "tools", "data", "sparse_linalg_heads.tsv"))
    ap.add_argument("--no-reference", action="store_true",
                    help="reuse --out as the reference table instead of re-running it")
    a = ap.parse_args()

    tpls = TEMPLATES
    if a.only:
        keep = set(a.only.split(","))
        tpls = [t for t in TEMPLATES if t[0] in keep]
    prog = program(tpls, SIZES, KINDS)
    tmp = tempfile.mkdtemp(prefix="sla_sem_")
    path = os.path.join(tmp, "probe.m")
    with open(path, "w") as f:
        f.write(prog)

    if a.no_reference:
        with open(a.out) as f:
            ref = parse(f.read())
    else:
        ws = find_wolframscript(a.wolframscript)
        if not ws:
            sys.exit("no wolframscript found")
        ref = parse(run([ws, "-file"], path, 3600))
        os.makedirs(os.path.dirname(a.out), exist_ok=True)
        with open(a.out, "w") as f:
            f.write("# id\tkind\tn\tsignature   (generated by tools/sparse_linalg_semantics.py)\n")
            for k in sorted(ref):
                f.write("\t".join(k) + "\t" + ref[k] + "\n")
        print(f"reference: {len(ref)} rows -> {a.out}")

    if a.mathilda:
        mine = parse(run([a.mathilda, "-file"], path, 3600))
        bad = 0
        for k in sorted(ref):
            if mine.get(k) != ref[k]:
                bad += 1
                print("\t".join(k), "ref=" + ref[k], "mathilda=" + str(mine.get(k)), sep="\t")
        print(f"{bad} of {len(ref)} rows differ")


if __name__ == "__main__":
    main()
