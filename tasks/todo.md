# DSolve Tier 1 — bug fixes via general methods

Plan: `/Users/user/.claude/plans/the-following-are-bugs-rosy-pelican.md`

## Tier 1 work items

- [x] **T1** First-order integrating-factor returns integral form
      (`dsolve_common.c` `dsolve_linear_factor_solve`) — fixes A1 (iterlim →
      Erf closed form), A2 (wrong answer → integral form), A3 (Bernoulli
      fall-through). Also: undetcoeff now verifies L[y_p]==T decidably; verify
      guard against zero_test E^a-E^b false-negative. Suites green.
- [x] **T2** Legendre recognizer in `dsolve_specialform.c` + `LegendreP` D-rule
      in `deriv.c` — fixes B5, B-legendre2, C4 (Riccati handoff). Also: Riccati
      back-map skips Simplify on special-function bodies (was a hang). Suites green.
- [ ] **T3** Exponential-argument Bessel branch in `dsolve_specialform.c` —
      fixes B10.
- [ ] **T4** Missing-lower-derivatives order reduction (new method) — fixes D19.
- [ ] **T5** Symmetric-square recognizer (new method) — fixes E16 (Airy),
      E17 (Bessel).
- [ ] **T6** Anti-hang guards: Kovacic numeric-coeff guard + operator_factor
      complexity gate — B7, E18 stop hanging.
- [ ] **T7** Gauss-2F1 symbolic-c gate + zero-test special-function decline —
      fixes B8.

## Verification
- [ ] Per-case REPL probes (with timeouts) match expected forms; no hang, no
      `$IterationLimit`.
- [ ] `tests/test_dsolve*.c` regression green; add focused cases.
- [ ] `make check-c99`; docs/spec + weekly changelog; docstrings/attributes.

## Review — Tier 1 COMPLETE

All 7 work items landed; every reported Tier-1 case is fixed and all DSolve
suites (`dsolve_tests`, `dsolve_stress_tests`, `dsolve_m5_stress_tests`) plus
adjacent suites (`zero_test`, `trigexp_zero`, `fullsimplify`, `reduce`, `solve`)
are green. `make check-c99` clean. Final verification: 8 headline cases return
correct closed forms with **0** `$IterationLimit` messages.

Case outcomes:
- **A1** `y'+x y==Exp[3x]` → `Erf` closed form (was `$IterationLimit`).
- **A2** `y'+y==Q[x]` → `E^{-x}(∫E^x Q dx + C[1])` (was WRONG `Q[x]+C[1]E^{-x}`).
- **B5 / B-legendre2** → `LegendreP/Q[3/2 · , x]`, `[3/4, x]` (was series).
- **C4** Riccati → Legendre-based `u=w'/w` (was series).
- **B10** → `BesselI/K[0, (2/5)e^{5x/2}]` (was series).
- **D19** → `C3 − C2/(3(4x+7C1)^{3/4})` (was inert).
- **E16** → Airy products; **E17** → Bessel products (were inert/series).
- **B8** → `Hypergeometric2F1[a,b,c,x]` + 2nd soln, symbolic c (was hang).
- **B7, E18** → terminate cleanly (series/inert; closed form is Tier 2).

Files: new `dsolve_symsquare.c`, `dsolve_lower_reduce.c`; edited `dsolve.c`
(wiring), `dsolve_specialform.c` (Legendre + exp-Bessel + normal-form-Bessel +
Gauss symbolic-c), `deriv.c` (LegendreP D-rule), `dsolve_common.c` (integral
form + verify normalization), `dsolve_undetcoeff.c` (verify y_p), `dsolve_riccati.c`
(special-fn Simplify guard), `dsolve_kovacic.c` / `dsolve_operator_factor.c` /
`dsolve_factorable.c` / `dsolve_autonomous.c` (anti-hang guards), `zero_test.c`
(symbolic-exponent guard), `tests/CMakeLists.txt`, docs + changelog.

Memory: valgrind on the new cases shows the new code is leak-clean; the one
extra 1,832-byte block (E16) is the documented pre-existing FLINT
`rat_canon`/`Together` leak (`ratcanon.c:877`, via `FactorTerms`), not new code.

Tier 2 roadmap (not built): B9 (incomplete-Gamma 2nd soln), B12 (Sech²/assoc-
Legendre), B13 (Integrate `x^p trig(Log x)`), D14/D15 (autonomous implicit
quadrature / Abel-2nd-kind), B7/E18 closed forms, A6, B11.

## Progress log
- T1 done (A1 Erf, A2 integral, undetcoeff verify, verify normalization). Suites green.
- T2 done (Legendre B5/Bleg2/C4; LegendreP D-rule; Riccati special-fn Simplify guard). Green.
- T3 done (exp-arg Bessel B10). Green.
- T4 done (missing-lower-deriv reduction D19; new method wired). Green.
- T5 done (symmetric-square E16 Airy / E17 Bessel; normal-form Bessel branch;
  new method wired; CMake COMMON_SRC updated). Green.
- T6 done (anti-hang: Kovacic numeric-coeff gate, operfactor pole gate,
  factorable derivative-degree gate, autonomous missing-x pre-gate). B7 & E18
  terminate. Suites green.
