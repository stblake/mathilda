# MonomialList / CoefficientRules / FromCoefficientRules

Plan: `/Users/user/.claude/plans/let-s-implement-coefficientrules-coeffic-fizzy-feather.md`
Sparse `{expvec -> coeff}` polynomial representation (+ monomial list + inverse).
Shared core in a new module; six named monomial orders + weight matrices; `Modulus`.

- [x] New `src/poly/monomials.c`: shared core (resolve vars, expand, term→expvec,
      merge, order) + `builtin_monomiallist` / `builtin_coefficientrules` /
      `builtin_fromcoefficientrules` + `monomials_init`.
- [x] Unified weight-matrix comparator (6 named orders synthesized + explicit matrix).
- [x] Wire: `poly.h` decls, `poly_init()` call, `sym_names.{h,c}` (3 symbols),
      `info.c` docstrings (3), `tests/CMakeLists.txt` COMMON_SRC + target.
- [x] Tests `tests/test_coefficient_rules.c` + CMake (all user examples, round-trips,
      order identities, argrx, Modulus, weight matrix). All green.
- [x] Docs: `structural-manipulation.md` (3 sections) + changelog `2026-09-14.md`.
- [x] Verify: build (clean), check-c99 (pass), suite (pass), valgrind (no leaks in
      new code — baseline is macOS libobjc/dyld), poly_tests (no regression),
      check-packed-aware (pass), check-nd-surfaces (pass), REPL (all examples exact),
      docstrings render via Information[].
- [x] check-fastpath-sweep flagged CoefficientRules/MonomialList (a packed vector
      handed as the vars arg is materialised to iterate as variables — inherent, no
      numeric surface); recorded both in nd_fastpath_sweep.py OFF_BUFFER with a reason.
      Re-run: my heads gone from the NEW list. The sweep's residual red is PRE-EXISTING
      (a rotating, nondeterministic backlog of unrelated heads — Exponent, SortBy, Div,
      Grad, ChineseRemainder, ... — from a stale OFF_BUFFER baseline last recorded
      2026-08-03), not caused by this change. Like check-compile-coverage's pre-existing red.

## Review

- Implemented all three heads in one self-contained module `src/poly/monomials.c`
  (~470 lines) sharing a single core; MonomialList/CoefficientRules differ only in
  output rendering, FromCoefficientRules is the inverse.
- Every user-supplied example reproduces Wolfram's output exactly, including the two
  algebraic order identities (2-var DegreeLex==DegreeRevLex; NegDegRevLex ↔
  reversed-vars DegreeLex) and both explicit weight-matrix examples.
- Design decision honoured: no-variable forms default to Mathilda's (sorted)
  `Variables[poly]`; the one adversarial case `CoefficientRules[y+x z]` uses `{x,y,z}`
  order (documented deviation from Wolfram's `{y,x,z}`).
- Symbolic coefficients handled via the Expr-based term walk (not the rational-only
  GBPoly path), so `a x + b x` merges to `{{1}->a+b}`.
- `Modulus->m` applied up front via the tested `PolynomialMod`; unknown order strings
  and bare List/NDArray inputs decline (return unevaluated); `FromCoefficientRules`
  wrong arg count emits `FromCoefficientRules::argrx` and stays unevaluated.
- Numeric-surface audits: symbolic/structural heads with no packed/NDArray/Compile
  surface by design; the static and dynamic ND audits do not flag them (no exempt-list
  entry needed). check-compile-coverage's red is pre-existing on main and does not name
  these heads.
- Gotcha caught by tests: CMake test build needs the new .c in COMMON_SRC (makefile
  auto-discovers; CMake does not). `CoefficientRules[poly, All]` correctly treats every
  Variables[poly] symbol (incl. a,b,c) as a variable — the initial test assumption was
  wrong, not the code.

---

# DSolve — M55 (§2.1.2 method wave) then M56 (M18 Stage 2)

Plan: `/Users/user/.claude/plans/what-is-the-next-shiny-bentley.md`
User direction: address option 2 (§2.1.2 method wave), then option 4 (M18 Stage 2),
with extensive unit + stress tests.

## M55 — generalised power-potential recogniser (§2.1.2 2nd-order-linear gap)

Refined target (the "generalized Bessel" plan framing was imprecise): the corpus
family `y'' + (α x^(2c) + β x^(c-1)) y == 0` is NOT plain Bessel — it reduces (via
ξ = x^(c+1)) to a **Coulomb/Whittaker** equation → `Hypergeometric1F1`. Only the
single-power `y'' + A x^m y == 0` case is Bessel. Both are now emitted, symbolic
exponent allowed, gated by `sf_num_ok` (0-FAIL by construction).

- [x] New `whittaker_M_1F1` + `specialform_power_potential` in `dsolve_specialform.c`
      (exponent grouping via `e = x T'/T`; `Expand` first so a factored normal-form
      potential splits into monomials; structural check `P_big − 2 P_sm == 2`).
- [x] Wire P==0 direct surface (self-gated by `sf_num_ok`).
- [x] Wire the Liouville normal-form pre-pass (covers `y'`-carrying members 803/804).
- [x] **Root-cause fix (dsolve_common.c `ds_residual_numeric_zero`)**: it declined the
      numeric KEEP whenever the residual had a `Derivative`/undefined head — but a
      `Derivative` of a DEFINED special function (Bessel/pFq/Airy) numericizes. That
      false decline routed symbolic-exponent special-function residuals to `zero_test`,
      whose precision ladder SPINS (name-sensitively — hung on exponent symbol `n`/`nn`,
      fine on `c`/`k`/`m`). Now declines only inert `Integrate`; arbitrary `f[x]` still
      yields NaN→skip→same zero_test path.
- [x] REPL spot-checks: 434, 792, 803, 804, single-power symbolic n all solve,
      residuals ~1e-15.
- [x] Unit test `t_m55_generalized_power_potential` (+ registered); all 4 assertions verified True.
- [x] Stress file `test_dsolve_m55_stress.c` (4 forward-generator families) + CMake; green.
- [x] m14/m18/m55 stress green (surgical fix has zero blast radius on non-x-power residuals).
- [x] §2.1.2 corpus re-measure: **564 PASS / 640 non-PASS / 0 FAIL** (+7 vs M54's 557); gate
      655 → 648; report/TSV regenerated; 9 flagship cases (72/81/103/434/791/792/803/804/805) PASS.
- [x] Docs: STATUS.md (M55 row + 2nd_linear 239/175), DSOLVE_PLAN.md M55, docs/spec/builtins/
      calculus.md + docstring + changelog; version 0.152 → 0.153.
- [x] check-c99 green; refresh code-review graph (incremental, all edits picked up).
- [~] valgrind spot-check (running — checking for leaks in specialform_power_potential/whittaker_M_1F1).

### Pre-existing unit-suite rot discovered (NOT M55; flagged for follow-up)
The `dsolve_tests` unit binary has been blocked at position 49 by the **documented pre-existing
rischnorman hang** (`t_rischnorman_enum_cap_no_crash`, SIGALRM 142) since ~M35, so tests after it
have not run to completion in days. Disabling that hang to run the full suite surfaced ≥2
**pre-existing** red tests (confirmed via A/B with the old verify — my M55 guard is a proven no-op
on their residuals, which carry no non-integer power of bare x):
- `t_kovacic_complex_poles`: `DSolve\`Kovacic[y''-((3+2x²)/(1+x²)²)y==0]` now DECLINES (a real
  Kovacic regression from some earlier milestone — masked).
- homogeneous IVP `y'==(x-y)/(x+y), y[1]==1`: now returns a CORRECT **implicit** solution
  (verified `{True}`) where the test expects an explicit form → stale test.
(There may be more after position 207.) All predate M55. The corpus gates (green, 0-FAIL) remain
the reliable functional check per project convention. **Recommend a separate cleanup pass** to
fix the rischnorman hang gating and refresh/repair the rotted post-49 unit tests.

## Kovacic ±i-pole fix (user-requested: "just the Kovacic regression") — DONE
- [x] Root-caused: `kovacic_case1_general` (`dsolve_kovacic.c`) blanket-declined ALL non-real
      poles (a guard vs a `ds_simplify(theta)` Heun hang, §2.2.14-1392/1393), which also broke the
      ±i case `y''−((3+2x²)/(1+x²)²)y==0` (masked for days by the pre-existing rischnorman unit hang).
- [x] Fix: decline a complex pole only when it has a NONZERO real part (√-discriminant radical =
      the spin); purely-imaginary pairs (`x²+c`) are solved. A/B-confirmed Heun cases still decline
      fast (no hang) → Frobenius.
- [x] `t_kovacic_complex_poles` passes (verified directly: `x√(1+x²)` etc., PossibleZeroQ True).
- [x] §2.1.2 re-measured: 564 → **568 PASS, 0 FAIL** (gate 648 → 644); §2.2.14 99/100, 0 FAIL,
      1392/1393 PASS via Frobenius. 2nd_linear 243/414 (58.7%).
- [x] Docs: changelog Kovacic entry, DSOLVE_PLAN M55 addendum, STATUS row; version 0.153 → 0.154;
      check-c99 green; graph refreshed.

## M56 — M18 Stage 2 (reducible-μ μ(x,y′)/μ(y,y′)) — implicit first-integral emission
(not started — begins after sign-off)

## Review
- **M55 DONE (bar valgrind confirm):** generalised power-potential recogniser (Bessel single-power
  + two-term Coulomb/Whittaker→1F1) at symbolic exponent, `sf_num_ok`-gated (0-FAIL). Root-cause
  verify fix (`ds_residual_numeric_zero` symbolic-x-power discriminator) removed a zero_test
  precision-ladder hang. §2.1.2 557→564 PASS, 0 FAIL, gate 655→648. v0.153. NOT committed.
