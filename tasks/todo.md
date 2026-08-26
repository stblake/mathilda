# Implement CharacteristicPolynomial

## Plan
Reuse eigen module internals. Ordinary case: Faddeev-LeVerrier (O(n^4)) + (-1)^n
sign fix. Generalized {m,a}: build (m - λa) + Laplace det. Substitute user var,
Expand. Returns symbolic polynomial (exempt from packed/Compile surfaces).

## Tasks
- [x] `src/linalg/charpoly.c` — new builtin `builtin_characteristicpolynomial`
- [x] `src/linalg/eigen.h` — declare builtin
- [x] `src/linalg/eigen.c` — register in `mateigen_init()` with ATTR_PROTECTED
- [x] `src/sym_names.h` / `src/sym_names.c` — add SYM_CharacteristicPolynomial (3 sites)
- [x] `src/info.c` — docstring in `info_init()`
- [x] `src/pack.c` — add to AWARE list (mirrors Eigenvalues; audit requires it)
- [x] `tests/CMakeLists.txt` — add source to COMMON_SRC + test executable block
- [x] `tests/test_characteristicpolynomial.c` — 19 test groups + leak loop
- [x] `docs/spec/builtins/linear-algebra.md` — entry
- [x] `docs/spec/changelog/2026-08-24.md` — changelog note
- [x] Build main binary clean, REPL smoke test (all reference cases match)
- [x] Build + run test suite (all pass), valgrind (no leak vs baseline)
- [x] `make check-c99`, `make check-packed-aware` — both clean

## Review

Implemented `CharacteristicPolynomial[m, x]` (`Det[m - x I]`) and the generalized
`CharacteristicPolynomial[{m, a}, x]` (`Det[m - x a]`) as a thin builtin over the
eigen module's existing char-poly machinery:

- **Ordinary case** → `eigen_char_poly_faddeev` (O(n^4)), negated by `(-1)^n`
  since it returns `det(λI - m)` but we want `Det[m - x I]`. This is what makes
  the 100×100 machine case fast (0.063s measured, ref ~0.09s) — the naive
  `Expand[Det[m - x I]]` would hit an O(n!) Laplace expansion of a symbolic-in-x
  matrix.
- **Generalized case** → `eigen_build_lambda_matrix` (m - λa) + `eigen_compute_det`
  (Laplace), correct sign directly. Shared null space → degree deficit (infinite
  generalized eigenvalue), verified (`-1 - x + x^2`, no x^3 term).
- Built in a private internal lambda, then the user's variable (symbol / number /
  expression) is substituted and the result Expand-ed.

All reference cases match (integer, symbolic, identity, zero, rational, machine,
complex, generalized, degree-drop, numeric-var, arity error). Returns a symbolic
`Plus` → exempt from packed/Compile surfaces but AWARE for NDArray *input*.

No corrections from the user during implementation → no new lessons.
