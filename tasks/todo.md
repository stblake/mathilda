# Implement `LegendreQ` (Legendre function of the second kind)

Plan: `/Users/user/.claude/plans/similar-to-our-implementation-valiant-haven.md`

Refinement vs plan: co-locate `LegendreQ` in `src/special_functions/legendre.c`
(reuses the static mpq/recurrence/poly helpers with zero duplication) instead of
a new `legendreq.c`. So no `core.c` / `COMMON_SRC` change; `legendre_init`
registers both.

## Tasks

- [x] `sym_names.h` / `sym_names.c` — declare/define/intern `SYM_LegendreQ`.
- [x] `legendre.c` — `builtin_legendre_q` + helpers (closed form, x=0 special
      value, non-integer numeric via Hypergeometric2F1 template, associated
      types 1/2/3), register in `legendre_init` with attributes.
- [x] `legendre.h` — declare `builtin_legendre_q`.
- [x] `sf_machine.c` / `.h` — `sf_machine_legendre_q` (recurrence, integer n, |x|<1).
- [x] `ndkernels.c` — `ndk_LegendreQ_c` + `NDKB_LegendreQ` + `REG_B`.
- [x] `deriv.c` — `SYM_LegendreQ` argument-derivative rule.
- [x] `info.c` — docstring.
- [x] `tests/test_legendreq.c` + `tests/CMakeLists.txt` (executable + add_test).
- [x] Docs: `docs/spec/builtins/special-functions.md` + `changelog/2026-08-31.md`.

## Verification
- [x] `legendreq_tests` all pass (15 groups; ctest-registered, unlike LegendreP).
- [x] REPL spot-checks match references (closed forms, 50-digit N, complex,
      associated, derivative, Series).
- [x] valgrind clean (byte-identical to startup baseline; zero new leaks).
- [x] Compile/NDArray parity (`Compiled -> True` scalar+array; packed==list).
- [x] `make check-c99`, `check-packed-aware`, `check-array-exactness` green.
- [x] `make check-nd-surfaces` green (exit 0; the 3 flagged heads
      minmax/extract/cross are pre-existing, not LegendreQ).
- [n/a] `make check-compile-coverage` fails on PRE-EXISTING image/array heads
      (ArrayPlot/Image*/PackedArrayQ, not in files I touched, not in BASELINE);
      LegendreQ itself compiles and is NOT in the failing list.

## Review

Implemented `LegendreQ` co-located in `legendre.c` (reuses the mpq recurrence /
`leg_poly_from_coeffs` / `leg_pm_power` machinery). Key results:
- Integer-order closed form `Q_n = P_n L + v_n` (v_n from the shared recurrence
  seeded 0, -1) — reproduces the reference `LegendreQ[0..5, x]`.
- Non-integer numeric via a `Hypergeometric2F1`/`Gamma` expression template —
  matches `N[LegendreQ[3/2,1/2],50]` to 50 digits, plus the complex case.
- Associated types 1/2/3 via `D[]` on `Q_n(t)` + branch-split prefactors —
  the type-2/3 construction reproduces `LegendreP`'s own outputs.
- `D[LegendreQ[n,x],x]` rule + origin `Series` (via naive Taylor + exact
  `Q_v(0)`).
- Machine kernel + binary ND kernel → packed/NDArray/Compile/auto-compile.
