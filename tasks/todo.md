# Plan: NIntegrate fixed-rule methods (Riemann / Trapezoidal / Newton–Cotes)

## Goal
Add three explicit `Method` strategies to `NIntegrate`, with method sub-options,
at both machine and arbitrary (MPFR) precision:

- `"RiemannRule"` with `"Type" -> "Left" | "Right" | "Midpoint"`
- `"TrapezoidalRule"` with `"RombergQuadrature" -> True | False` (default `True`)
- `"NewtonCotesRule"` with `"Points" -> n` (default 3 = Simpson)

Today `NI_TRAP` exists in the enum but is unimplemented (warns + unevaluated);
`"RiemannRule"`/`"NewtonCotesRule"` are unrecognized (`NI_UNIMPL`). Method
sub-options are currently discarded.

## Design

All three are equally-spaced composite rules sharing one refinement loop:
build a base estimate, halve the step, optionally Richardson-extrapolate
(Romberg), stop when the successive-estimate error meets PrecisionGoal/AccuracyGoal
or a level/eval cap is hit. Endpoints **are** sampled (matches WL for these rules).
Accumulate in `double _Complex` (machine) / `mpfr_t` re-im pairs (MPFR), matching
the `denint` convention.

- **Riemann** Left/Right/Midpoint: first-order; refine by doubling panels,
  error = |Iₙ − I₂ₙ| (no Romberg). Reuse prior samples where possible.
- **Trapezoidal**: standard recursive trapezoid (reuses midpoints). Romberg=True
  → Richardson table T[j][m]; best = T[j][j]. Romberg=False → best = T[j][0].
- **Newton–Cotes**: composite closed N–C of order `Points` (2=trap, 3=Simpson,
  4=Simpson-3/8, 5=Boole); refine by doubling panels + Romberg extrapolation.

## Steps

### 1. New kernel: `src/numerical_calculus/ncrule.{h,c}`
- `ncrule.h`: include `gkadapt.h` (reuse `GkSampleMachine`); enum
  `NcrRuleKind { NCR_RIEMANN_LEFT, _RIGHT, _MIDPOINT, NCR_TRAPEZOIDAL, NCR_NEWTONCOTES }`;
  declare `ncr_integrate_machine(...)` and (USE_MPFR) `ncr_integrate_mpfr(...)`.
- `ncrule.c`: implement both. Args: callback+ctx, a, b, kind, nc_points, romberg,
  reltol, abstol, max_levels, out value + abserr (machine: `double _Complex*`;
  MPFR: `mpfr_t` re/im + bits, mirroring `denint_tanhsinh_mpfr`). Strict C99,
  guarded `M_*` if any, careful `mpfr_init2`/`mpfr_clears` pairing.

### 2. Wire into `src/numerical_calculus/nint.c`
- `#include "ncrule.h"`.
- `NiMethod`: add `NI_RIEMANN`, `NI_NEWTONCOTES` (keep `NI_TRAP`).
- `ni_method_implemented()`: add the three.
- `ni_method_from_string()`: add `"RiemannRule"`→`NI_RIEMANN`,
  `"NewtonCotesRule"`→`NI_NEWTONCOTES` (TrapezoidalRule already → `NI_TRAP`).
- `NiOpts`: add `int rule_type` (0=Left default), `bool romberg` (default true),
  `int nc_points` (default 3). Set defaults at init in `builtin_nintegrate`.
- `ni_apply_option()` Method branch: iterate remaining `List` args; parse
  `"Type"->str`, `"RombergQuadrature"->True/False`, `"Points"->int` (string LHS,
  string/bool/int RHS). Unknown sub-opts: ignore (or warn `badopt`).
- `ni_core_finite()`: add branches calling a new `ni_try_ncr(...)` for
  `NI_RIEMANN`/`NI_TRAP`/`NI_NEWTONCOTES`, mapping `rule_type`→`NcrRuleKind` for
  Riemann. Cap iterations via a level/MaxPoints budget (Riemann needs more).
- MPFR path in `ni_run_1d_finite_real()`: when method is a fixed rule, call
  `ncr_integrate_mpfr` instead of `denint_tanhsinh_mpfr`; else unchanged.

### 3. Docstring
- Add a terse `symtab_set_docstring("NIntegrate", ...)` in `nintegrate_init()`
  (none exists today — pre-existing gap). No examples (per project rule).

### 4. Build wiring
- `makefile`: auto-discovered (no change).
- `tests/CMakeLists.txt`: add `../src/numerical_calculus/ncrule.c` to COMMON_SRC
  (~line 333, beside `nint.c`) so `nint_tests` links.

### 5. Tests — `tests/test_nint.c`
Add `ASSERT_CLOSE` cases (the user's examples) + Midpoint + an MPFR case:
- Riemann Left @PrecisionGoal 2 → 12.1732; Right → 12.0847
- Trapezoidal (RombergQuadrature->False) → 12.1561; default (Romberg) → 12.1561
- NewtonCotesRule → 12.1561
- `WorkingPrecision->30` TrapezoidalRule on a smooth integral
  (`Exp[Cos[x]]`,{x,0,10}) vs the DoubleExponential value, to exercise MPFR.

### 6. Docs + changelog
- `docs/spec/builtins/numerical-calculus.md`: add the three to the named-methods
  list; document the sub-options; remove `"NewtonCotesRule"` from the
  not-implemented example.
- `docs/spec/changelog/2026-06-15.md` (Monday of current ISO week): add a summary.

### 7. Verify
- `make -j` clean build (foreground; C99 `-Wall -Wextra`).
- Build + run `nint_tests` (foreground, scoped to this binary).
- Manual REPL check of all six example invocations + the MPFR case.
- `valgrind` spot-check on a couple of the new method calls (diff vs baseline
  noise per project memory).

## Review

Implemented and verified.

- **New kernel** `src/numerical_calculus/ncrule.{c,h}`: machine + MPFR composite
  Riemann (Left/Right/Midpoint), Trapezoidal (Romberg on/off) and Newton–Cotes
  (Points 2–5) rules. Two-row Romberg (only level j−1 kept) bounds MPFR alloc.
- **`nint.c`**: `NI_RIEMANN`/`NI_TRAP`/`NI_NEWTONCOTES`, Method sub-option
  parsing (`"Type"`/`"RombergQuadrature"`/`"Points"`), `ni_try_ncr` dispatch in
  `ni_core_finite`, MPFR-path routing, plus the missing `NIntegrate` docstring.
- **Key bug found & fixed**: the Newton–Cotes *MPFR* path multiplied samples by
  double-precision weights (4/3, 1/45, …), capping accuracy at ~1e-16 regardless
  of `WorkingPrecision`. Replaced with exact integer numerators + one final
  division by the shared denominator → Simpson now exact for cubics at full
  precision; Trapezoidal gives Pi to 38 digits at WP40.
- **Verification**: `make` clean, `-Wall -Wextra` zero warnings. `nint_tests`
  all 24 groups pass (added 4 new test groups, ~45 new assertions). valgrind
  clean (lost = documented macOS baseline 12,800B/400 blocks; no ncrule/nint
  frames in any lost stack). Docs + changelog updated.
- **WL fidelity caveat**: convergent methods (Trapezoidal/Newton–Cotes) match WL
  exactly (12.1561). WL's crude-precision *Riemann* values come from plugging the
  rule into its adaptive (non-uniform) subdivision strategy; we implement the
  textbook uniform composite rule, so crude Riemann values differ (documented).
