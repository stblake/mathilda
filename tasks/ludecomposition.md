# LUDecomposition Implementation Plan

## Goal

Implement `LUDecomposition[m]` returning `{lu, p, c}` where:
- `lu` is the combined L/U matrix (Doolittle form: L unit lower
  triangular, U upper triangular, stored compactly in one n×n matrix);
- `p` is a 1-indexed row-permutation vector with `m[[p]] == l . u`
  (where `l = LowerTriangularize[lu, -1] + IdentityMatrix[n]` and
  `u = UpperTriangularize[lu]`);
- `c` is an L∞-condition-number estimate for approximate numerical
  matrices, or exact Integer `0` for exact / symbolic inputs.

Follows the QRDecomposition split: top-level dispatcher routes by
leaf precision to a symbolic core (Doolittle through the evaluator),
a LAPACK machine-precision kernel, or an MPFR kernel; on any failure
the upper layers fall back to the symbolic core.

## Files

- `src/sym_names.{c,h}` — add `SYM_LUDecomposition`.
- `src/linalg/ludecomp.h` — public `builtin_ludecomposition` and `ludecomp_init`.
- `src/linalg/ludecomp_internal.h` — shared kernel signatures.
- `src/linalg/ludecomp.c` — builtin entry, dispatcher, symbolic core.
- `src/linalg/ludecomp_machine.c` — LAPACK kernel (`dgetrf`/`zgetrf` +
  `dgecon`/`zgecon` + `dlange`/`zlange`).
- `src/linalg/ludecomp_mpfr.c` — MPFR Doolittle kernel + L∞ condition
  number from the explicit inverse.
- `src/linalg/lapack.{h,c}` — add Fortran prototypes and thin C
  wrappers for `dgetrf`, `zgetrf`, `dgecon`, `zgecon`, `dlange`,
  `zlange` (plus USE_LAPACK=0 stubs).
- `src/core.c` — call `ludecomp_init()`.
- `src/info.c` — `LUDecomposition` docstring.
- `docs/spec/builtins/linear-algebra.md` — `## LUDecomposition` section.
- `docs/spec/changelog/2026-05.md` — changelog entry.
- `tests/CMakeLists.txt` — register three test binaries; add the new
  source files to `COMMON_SRC`.
- `tests/test_ludecomposition.c` — symbolic + exact integer / rational
  / complex.
- `tests/test_ludecomposition_machine.c` — random doubles + condition
  number sanity.
- `tests/test_ludecomposition_mpfr.c` — high-precision MPFR pieces.

## Algorithm (Doolittle with partial pivoting)

```
LU = A  (in-place); perm = [1..n]; sign = +1
for k = 0 .. n-1:
    pivot_row = k
    if numeric:                # largest |A[i,k]|
        pick the row in [k, n) maximising |A[i, k]|
    else (symbolic / exact):    # first non-zero pivot from row k down
        if A[k,k] is definitely zero:
            scan rows below for first non-zero
    swap rows k and pivot_row in LU, swap perm[k]/perm[pivot_row]
    pivot = LU[k, k]
    if pivot is non-zero:
        for i = k+1..n-1:
            LU[i, k] /= pivot
            for j = k+1..n-1:
                LU[i, j] -= LU[i, k] * LU[k, j]
```

Symbolic kernel runs every primitive through the Mathilda evaluator
(`eval_plus`, `eval_times`, `eval_power`, `Together` for the
zero-detection guard, `Expand` for tidy-up); inexact inputs use the
`PseudoInverse` rationalise → exact → numericalise round-trip so the
symbolic core also covers Real / MPFR matrices when the fast kernels
aren't selected.

## Condition number

- Exact / symbolic input: emit `Integer 0` (matches the Mathematica
  example where exact LU returns `c = 0`).
- Machine-precision: use `dgecon` (`'I'` norm) with `anorm = dlange('I')`,
  return `1.0 / rcond` as `EXPR_REAL`.
- MPFR: invert L and U via back-substitution, multiply to get A^{-1},
  return `||A||_∞ * ||A^{-1}||_∞` as `EXPR_MPFR` at the input
  precision.  This is O(n^3) which matches the cost of the
  factorisation itself.

## Singular handling

When the kernel discovers a zero pivot it emits
`LUDecomposition::sing: Matrix m is singular.` (one-shot via
`matsol_warn_once`) and *continues* the factorisation, leaving the
zero on the diagonal of U — matching Mathematica's behaviour shown
in the prompt's example.

## Pivoting convention

`p[k]` is the 1-indexed original-row number that was used as the
pivot at step `k`.  Equivalent to LAPACK's `ipiv` interpreted as
"the inverse of the row permutation: A_perm = A[p, :]".

## Memory contract

Standard builtin contract per SPEC.md §4: the builtin entry returns
either a freshly-allocated `Expr*` (success) or `NULL` (cannot
evaluate, evaluator owns `res`).  Every kernel uses the
NULL-out-before-free pattern for any reused sub-trees.  Each test
binary is intended to run cleanly under
`valgrind --leak-check=full`.

## Phases (= TaskCreate items)

1. Add `SYM_LUDecomposition` + register builtin.
2. Add LAPACK wrappers (`getrf`, `gecon`, `lange`).
3. Implement `ludecomp.c` (symbolic core + dispatcher).
4. Implement `ludecomp_machine.c` (LAPACK kernel).
5. Implement `ludecomp_mpfr.c` (MPFR kernel).
6. Add docstring + spec docs + changelog.
7. Wire build + write unit tests.
8. Build, run tests, valgrind smoke.

## Review (filled after completion)

All eight phases are complete and all three test binaries pass:

- `ludecomposition_tests` -- 14 symbolic / exact / complex / singular
  tests, including a print-form regression for `LUDecomposition[{{a, b},
  {c, d}}]` and an exact match for the Mathematica documentation
  Vandermonde example `LUDecomposition[{{1, 1, 1}, {2, 4, 8}, {3, 9,
  27}}]`.
- `ludecomposition_machine_tests` -- 8 LAPACK-path tests covering
  3x3 / 4x4 / 5x5 real, 2x2 / 3x3 complex, identity, singular, and
  mixed Integer/Real inputs.  The condition-number estimate for the
  Mathematica doc example `{{1.6, 2.7, 3.6}, {1.2, 3.2, 5.2}, {3.3,
  3.4, 6.5}}` comes out at exactly `20.8391`, matching the doc.
- `ludecomposition_mpfr_tests` -- 6 high-precision tests at 60, 80,
  100, and 200 bits.  Residuals scale with precision (`1.7e-80` at
  80 bits, `6.5e-201` at 200 bits) confirming the MPFR kernel is
  genuinely using the requested precision.

`valgrind --leak-check=summary` shows the same "definitely lost"
baseline as a no-LU run (`12,800 bytes in 400 blocks`), i.e. the new
code adds no leaks.

QR, linalg, nullspace, and matrank suites still pass; no regressions
from the LAPACK wrapper additions.

Two findings worth recording in `tasks/lessons.md`:

1. Mathilda's `Module` does NOT substitute its local variables into
   `HoldAll` bodies (`Table` in particular), so a Module-bound `lu`
   is left as the literal symbol inside `Table[lu[[i, j]], ...]`.
   `Block` works correctly for this pattern.

2. Picking iteration-variable names that collide with matrix entries
   (e.g. `Table[..., {i, n}, {j, n}]` while the matrix contains a
   symbol called `i`) silently corrupts the output -- `i` from the
   matrix takes precedence in evaluation. Use `ii` / `jj` (or any
   unambiguous name) for inner-loop iteration variables in tests
   that exercise general symbolic matrices.
