# Vector Analysis: Grad, Div, Curl, Laplacian

## Plan
Add four Protected builtins in `src/vectoranal.c`, built on top of `D`.
Cartesian forms (fully general rank) + chart forms (Cartesian/Polar/Cylindrical/Spherical).

## Tasks
- [x] Add SYM_Grad/Div/Curl/Laplacian to sym_names.h + sym_names.c
- [x] Create src/vectoranal.h
- [x] Create src/vectoranal.c (Grad, Div, Curl, Laplacian; Cartesian + charts)
- [x] Wire vectoranal_init into core.c core_init
- [x] Add docstrings in info.c
- [x] Build main binary, REPL smoke test against user examples (24/24 correct)
- [x] Create tests/test_vectoranal.c
- [x] Register test in tests/CMakeLists.txt (COMMON_SRC + 3-line target)
- [x] Run unit tests (ctest) — vectoranal_tests + deriv/deriv_array regression PASS
- [x] valgrind leak check — no leaks beyond core_init baseline
- [x] make check-c99 — PASS
- [x] Docs: docs/spec/builtins/calculus.md + changelog 2026-08-17

## Review

Implemented four Protected builtins in `src/vectoranal.c`, all assembled from
`D[...]` and reduced with one `evaluate()` (no new differentiation code):

- **Grad** = `D[f, {{vars}}]` (scalar→vector, vector→Jacobian, tensor→+1 rank).
- **Div** contracts the innermost slot (vector→scalar, tensor→rank-1 map).
- **Curl** = generalized Levi-Civita permutation sum (2-D scalar, 3-D vector,
  rank-2 tensor→scalar all one path).
- **Laplacian** = `Sum_i D[f,{x_i,2}]` (element-wise over arrays).
- **Charts**: Cartesian/Polar/Cylindrical/Spherical via scale factors
  (scalar Grad, vector Div, scalar Laplacian, 2-D/3-D vector Curl). Unsupported
  ranks + unknown charts left unevaluated (unknown chart emits `Head::chart`).

Verified against Wolfram reference outputs: Cartesian forms match exactly
(incl. rank-2 curl `1/2(x^3-3xy^2-3x^2y^2+2xy^3)`); chart forms reproduce all
Polar/Cylindrical/Spherical scalar & vector examples. Zero new memory leaks
(baseline-identical valgrind), C99-clean, no regressions in D suites.

Files: new `src/vectoranal.{c,h}`, `tests/test_vectoranal.c`; edited
`src/sym_names.{h,c}`, `src/core.c`, `src/info.c`, `tests/CMakeLists.txt`,
`docs/spec/builtins/calculus.md`, `docs/spec/changelog/2026-08-17.md`.
