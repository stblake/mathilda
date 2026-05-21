# QRDecomposition Implementation Plan

## Goal

Implement `QRDecomposition[m]` returning `{q, r}` such that
`m == ConjugateTranspose[q] . r` with `q` row-orthonormal (or
row-unitary in the complex case) and `r` upper-triangular ("thin"
form: `Length[q] == Length[r] == MatrixRank[m]`).

Support:
- exact symbolic matrices (integer / rational / Complex / free
  symbols)
- machine-precision real matrices (Real entries) - exact-pipeline
  with rationalize / numericalize round-trip so we get inexact-in /
  inexact-out at machine precision.
- arbitrary-precision (MPFR) matrices - same exact pipeline at the
  input bit-precision.
- non-square (n > p, p > n) and rank-deficient inputs.
- `Pivoting -> True` returning `{q, r, p}` with `m . p == q^T . r`.
- `TargetStructure -> "Dense"` (default) - emits plain lists. The
  `"Structured"` value is left unevaluated for now (it would require
  new `OrthogonalMatrix` / `UnitaryMatrix` / `UpperTriangularMatrix`
  / `PermutationMatrix` builtin heads we do not yet have).

## Algorithm

Compute the standard "thin" QR

    A = Q R

with Q (n x rank) having orthonormal columns and R (rank x p) upper-
trapezoidal, then return `q = ConjugateTranspose[Q]` and `r = R`.
This satisfies `m == ConjugateTranspose[q] . r` because
`ConjugateTranspose[ConjugateTranspose[Q]] == Q`.

Single implementation: **Modified Gram-Schmidt** driven through the
Mathilda evaluator so symbolic, exact rational, and complex paths
share one code path. For inexact (Real / MPFR) inputs we follow the
`PseudoInverse` pattern - `common_rationalize_input` first, run the
exact pipeline, then `common_numericalize_result` at the input
precision. This delivers machine-precision and MPFR semantics with
no separate Householder kernel and zero precision-floor surprises.

### Pivoting

When `Pivoting -> True` we maintain a permutation `perm` and at each
step choose the pivot column with the largest residual squared-norm
(after projecting away previously-built orthonormal vectors). For
exact inputs this still works - we use `expr_compare` on the
squared norms after `Together` / `Expand`. The resulting permutation
matrix `p` is built from `perm` so that column `j` of `m.p` is
column `perm[j]` of `m`, and `m . p == ConjugateTranspose[q] . r`.

### Rank-deficient handling

When the orthogonalised residual vanishes (`is_zero_poly` of the
squared norm), the column is in the span of the orthonormal vectors
built so far; we record its `r[:, k]` coefficients but do not add a
new row to `q` or `r` ("thin" QR).

## File layout

- `src/linalg/qrdecomp.h` - public surface (`builtin_qrdecomposition`,
  `qrdecomp_init`).
- `src/linalg/qrdecomp.c` - implementation.
- `src/sym_names.{c,h}` - intern `QRDecomposition`, `Pivoting`,
  `TargetStructure`.
- `src/core.c` - call `qrdecomp_init()` from the linear-algebra
  block.
- `src/info.c` - docstring (`Information[QRDecomposition]`).
- `tests/CMakeLists.txt` - add `qrdecomposition_tests` target.
- `tests/test_qrdecomposition.c` - unit tests.
- `docs/spec/builtins/linear-algebra.md` - public docs entry.
- `docs/spec/changelog/2026-05.md` - summary entry.

## Memory contract

Standard builtin contract (per `MEMORY.md` / SPEC.md §4.1): the
builtin owns `res` on success path implicitly (the evaluator frees
it on a non-NULL return), and the builtin must NOT itself call
`expr_free(res)`. All intermediate `eval_and_free` calls consume
their argument trees; we explicitly free every other allocation.

## Tests

- 2x2 invertible exact integer matrix - check `ConjugateTranspose[q].r == m`.
- 3x2 exact integer - check shape (q is 2x3, r is 2x2) and product.
- 2x3 exact integer - check shape (q is 2x2, r is 2x3) and product.
- 3x3 rank-2 (the famous `{{1,2,3},{4,5,6},{7,8,9}}`) - q and r are 2x3.
- 5x5 random-like exact - q.q^T == I_5 after Simplify; r upper triangular.
- Symbolic 2x2 `{{a,b},{c,d}}` - q^T.r equals m.
- Complex 3x3 - q.q^H == I; m == q^H.r.
- Real machine-precision 3x3 - check norms via Chop.
- MPFR matrix at 50 bits - check the precision flows through.
- `Pivoting -> True` on a 4x4 - first diagonal of r has largest magnitude.
- Edge cases: zero matrix returns `{{}, {}}`; 1x1 returns trivial QR.
- valgrind: every test run leak-free under `--leak-check=full`.
