# JordanDecomposition[m]

Return `{s, j}` with `m == s . j . Inverse[s]`; `j` = Jordan canonical form.

## Plan
- [ ] `src/linalg/jordandecomp.{c,h}` — builtin + exact/symbolic chain engine + numeric fast path
- [ ] Register in `core.c` (`jordandecomp_init`), attribute `Protected`
- [ ] Docstring in `src/info.c` (no examples)
- [ ] `src/pack.c` AWARE list += "JordanDecomposition" (not INT64_OK)
- [ ] `tools/compile_coverage.py` EXEMPT += "JordanDecomposition" ("returns a pair of matrices")
- [ ] `docs/spec/builtins/linear-algebra.md` — `## JordanDecomposition`
- [ ] `docs/spec/changelog/2026-08-24.md` — changelog entry
- [ ] `tests/test_jordandecomp.c` + `tests/CMakeLists.txt` (COMMON_SRC + test target)
- [ ] Build, run tests, valgrind, `make check-*`

## Algorithm
- Exact/symbolic: charpoly (`eigen_char_poly_faddeev`) → `eigen_solve_poly` → distinct
  eigenvalues with multiplicity. Per λ: nullity sequence of (m−λI)^k via
  `eigen_null_space`; top-down chain-top selection (extend-to-basis via MatrixRank),
  bottom-up chains. Assemble S (chains as columns), J (blocks). Gate: total cols == n
  else NULL (irrational-defective limitation).
- Numeric (inexact): fast path = numeric `Eigenvectors` as columns of S + diagonal J
  from component-ratio eigenvalues, when MatrixRank[evecs]==n (diagonalizable).
  Else rationalize → exact core → numericalize.

## Plan (done)
- [x] `src/linalg/jordandecomp.{c,h}` — builtin + exact/symbolic chain engine + numeric fast path
- [x] Register in `core.c` (`jordandecomp_init`), attribute `Protected`
- [x] Docstring in `src/info.c` (no examples)
- [x] `src/pack.c` AWARE list += "JordanDecomposition"
- [x] `tools/compile_coverage.py` EXEMPT += "JordanDecomposition"
- [x] `docs/spec/builtins/linear-algebra.md` — `## JordanDecomposition`
- [x] `docs/spec/changelog/2026-08-24.md` — entry
- [x] `tests/test_jordandecomp.c` (20 checks) + CMake target
- [x] Build clean, tests pass, valgrind clean, audits

## Review

**Result.** `JordanDecomposition[m]` → `{s, j}`, `m == s.j.Inverse[s]`. All spec
examples reproduced: exact `j` for the two 3×3s (`{{6,0,0},{0,12,1},{0,0,12}}`,
`{{24,0,0},{0,48,1},{0,0,48}}`), the size-3 chain 4×4, the 2×2 symbolic
(`m.s==s.j`), machine real/complex, MPFR (20-digit), the defective-`N` block,
and 1×1. 20/20 unit tests pass; eigen/linalg/matinv/nullspace suites unregressed.

**Design.** One field-agnostic exact chain engine (charpoly → distinct
eigenvalues with multiplicity → per-λ nullity sequence of `(m−λI)^k` via
`eigen_null_space` → top-down chain-top selection with a `MatrixRank`-based
extend-to-basis → bottom-up chains). Numeric fast path: distinct spectrum ⇒
diagonalizable ⇒ eigenvectors-as-columns + `DiagonalMatrix[eigenvalues]` (paired
positionally — verified reliable to ~1e-14), no inverse/product formed, so
100×100 runs in ~15 ms. Repeated numeric eigenvalue ⇒ rationalize → exact core →
numericalize.

**Perf note.** The numeric path deliberately avoids `Inverse`/`Dot` on the
(generally complex) eigenvector matrix: Mathilda has no packed complex linear
algebra, so boxed complex `Inverse[100×100]` is ~15 s. The distinctness gate
sidesteps it entirely.

**Surfaces.** AWARE (not INT64_OK); `Compile[]` exemption documented. c99 +
packed-aware audits green; `JordanDecomposition` not flagged by
compile-coverage (its pre-existing 22-head Image*/Interpolation backlog is
unrelated and not a CI gate). Valgrind: 160 calls add 0 lost bytes over baseline.

**Limitation (documented).** Exact matrix with an *irrational, defective*
eigenvalue is left unevaluated (the `is_zero_poly` pivot test can't span the
generalized eigenspace) rather than returned wrong.
