# P1 — Structure-theorem decision + logarithmic-derivative decision

Roadmap item P1 from `RISCH_BRONSTEIN_GAPS.md`. Bronstein §9.3 (Cor. 9.3.1/9.3.2),
§5.12 (In-Field Integration), §7.3 (Parametric Log Derivative). Substrate
(`risch_structure.c` complex case + `risch_rational_span`) already landed.
User-approved scope: **A + B + C-full (full oracle replacement)**. Extensive unit +
stress tests for each Bronstein routine. (Book page N ≈ PDF page N+17.)

## Phase A — Real structure case (Cor. 9.3.2 iii/iv)  [risch_structure.c]  ✅ DONE
- [x] Add `ArcTan` monomial kind to `rs_decode_tower`; expose per-monomial kind.
- [x] Partition generators by the four DISJOINT index sets E/L/T/A (E∪L for
      log/exp, T∪A for tan/arctan). Map coefficients back to full monomial order.
- [x] `Risch`TanReducible` (eq. 9.15) + `Risch`ArcTanReducible` (eq. 9.14).
- [x] Docstrings + attributes + changelog.
- [x] `tests/test_risch_structure_real.c` — all green (real reductions, genuine-new
      declines, disjoint-index, deeper real towers, back-compat, robustness).

## Phase B — LogarithmicDerivativeOfRadical (§5.12 / eq. 7.44)  [risch_structure.c]  ✅ DONE
- [x] `Risch`LogarithmicDerivativeOfRadical[f, x, mons]` -> `{n, u}` or `False`,
      structure-theorem test (Cor. 9.3.1/9.3.2 ii) via E∪L span + witness
      u = prod base^{n r}, D[u]/u == n f. (Full §5.12 recursion w/ log-monomial
      radicals = documented completeness item, pinned as sound decline.)
- [x] Docstrings/attrs/changelog.
- [x] `tests/test_risch_logderiv_radical.c` — all green (structural witnesses,
      end-to-end D[u]/u−n f==0 identity, declines, §5.12 scope boundary, stress).

BASELINE (pre-Phase-C, must not regress): integrate_risch_transcendental_tests &
integrals_tests CLEAN; intrat_tests 2 FAIL + intrischnorman_tests 6 FAIL are
PRE-EXISTING (files untouched).

## Phase C — Full oracle replacement in the live integrator  [ATTEMPTED, REVERTED]
Outcome (2026-07-14): NET-NEGATIVE, reverted (integrator back to baseline).
- Blanket structure-oracle log collapse (rt_log_oracle_normalize pre-pass)
  REGRESSED the Bronstein `Log[x/Log[x]]` example: a composite log inside a
  `1/Log[...]^2` denominator becomes a COUPLED two-log denominator the field
  integrator declines. Collapsing only helps for redundant EXTRA generators, not
  kernels in denominators.
- Probe confirmed the live engine already handles dependent generators robustly
  (each becomes its own tower variable + triangular Dt, gated by diff-back); the
  elementary additive-dependent cases already integrate. G-A2's live impact is
  overstated. See memory project_risch_tower_dependent_log_robust.
- Decision: P1 ships as A+B (the reusable Bronstein decision builtins). A safe
  structure-oracle integration would need denominator-aware basis selection
  INSIDE rt_tower_build — a larger, delicate effort with unclear marginal benefit
  given diff-back already delivers soundness. Deferred as a future surgical item.

## Verification (each phase)
- Clean `-std=c99 -Wall -Wextra`; scoped `*_tests` green (grep FAIL:); REPL pins;
  valgrind vs Sin[1.0] baseline. Commit per phase.

## Cherry substrate refactors (CHERRY_DESIGN.md §3) — behaviour-preserving  ✅ DONE
- [x] R5: `PolynomialSqrt[p]`/`[p,x]` in facpoly.c + test_polynomialsqrt.c (green, leak-clean).
- [x] R2: `RtSpecialForm` registry for `rt_special_case` — byte-identical (regression clean).
- [x] R3: `Integrate`RothsteinTragerResultant[num,den,z,x]` in intrat.c (thin over `Resultant`)
      + test_rt_resultant.c (green, leak-clean).
- [x] R1: `rt_tower_solve` extension-point contract documented (already general — no change).
- [x] R4: DEFERRED to the Cherry module (internal recursion entry `rt_field_integrate` already
      exists; a wrapper now = dead code). Documented in CHERRY_DESIGN.md §3.

## Review
P1 landed on main (real structure theorem + LogarithmicDerivativeOfRadical); the wholesale
live-oracle replacement (P1 Phase C) was attempted, proven net-negative (regresses Bronstein
Log[x/Log[x]]), and reverted. Cherry incorporation was DESIGNED (CHERRY_DESIGN.md, grounded in
both 1986/1989 papers) and its 5 substrate refactors landed: 2 new tested primitives
(PolynomialSqrt, RothsteinTragerResultant), 1 registry (RtSpecialForm), 1 contract doc
(rt_tower_solve), 1 principled defer (lower-field hook). No current integral changed;
pre-existing failures (intrat 2 / factor 2 / poly 2 / intrischnorman 6) unchanged.
