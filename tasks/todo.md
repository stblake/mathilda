# Task: Implement LegendreP

## Scope (faithful, verifiable core)
- [ ] `src/special_functions/legendre.{c,h}` — new module
- [ ] `LegendreP[n, x]`:
  - exact integer n (any sign, P_{-1-n}=P_n) -> exact monomial polynomial (GMP rationals)
  - `LegendreP[n, 1] = 1` (any n)
  - non-integer n with an inexact arg -> numeric Gauss 2F1(-n,n+1;1;(1-x)/2), real+complex, arbitrary precision
  - else symbolic (NULL)
- [ ] `LegendreP[n, m, x]` (type 1, default): integer n,m>=0 -> (-1)^m (1-x^2)^(m/2) D^m P_n(x); m>n -> 0
- [ ] `LegendreP[n, m, a, x]`: a in {1,2,3}; a=1 -> type 1; a=2/3 -> 2F1Reg core * prefactor
- [ ] Attributes: Listable, NumericFunction, Protected
- [ ] sym_names.{c,h}: SYM_LegendreP
- [ ] core.c: include + legendre_init()
- [ ] info.c: docstring
- [ ] tests/CMakeLists.txt: COMMON_SRC + legendre_tests target
- [ ] tests/test_legendre.c: extensive
- [ ] docs/spec/builtins + changelog
- [ ] valgrind clean

## Deferred (documented): symbolic Series/SeriesCoefficient, D rules, non-integer associated/type forms, |w|>=1 analytic continuation.

## Review
DONE. All checkboxes complete.
- legendre.{c,h} implemented; exact integer polynomials via GMP-rational
  three-term recurrence; numeric non-integer order via Gauss 2F1 on ncpx
  (real/complex, machine + MPFR, verified to 300 digits vs WL); associated
  type-1 via Rodrigues derivative; types 2/3 via regularized-2F1 core.
- All WL example outputs reproduced exactly (poly forms, 9.58312, complex
  5.20466+0.299479 I, N[..,50], 300-digit timing, 10/2/x, types).
- Attributes {Listable, NumericFunction, Protected}; docstring; SYM_LegendreP.
- tests/test_legendre.c (13 groups) passes, exit 0, no FAIL.
- valgrind: identical to Sin[1.0] baseline (12,800B/400 — known macOS noise);
  zero LegendreP frames.
- Strict C99 -Wall -Wextra clean.
- Deferred (left symbolic, documented): non-integer associated/type forms,
  negative m, Series/SeriesCoefficient, D[] rules, |(1-x)/2|>=1 continuation.
