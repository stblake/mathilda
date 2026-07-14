# Risch Transcendental → Complete Bronstein Decision Procedure

Plan: `RISCH_COMPLETENESS_PLAN.md`. Reference: Bronstein, *Symbolic Integration I*, 2nd ed.
Each item lands green (its tests + the full `integrate_risch_transcendental_tests` suite) and
is committed before the next. Correct-by-construction: emit only behind an exact tower-var
`Together`-zero identity; every NULL is authoritative or a sibling-decline.

## Gap 1 — Recursive Risch DE over the tower (KEYSTONE)

Replace the `rt_field_rde` general-branch ansatz (L4122–4280) with a literal recursive Ch.6
solver `rde_tower(f, g, RdeCtx*)` solving `D_tower[y] + f y = g` over `K_L`.

- [x] **1a** — `RdeCtx` derivation abstraction + `rde_tower` + base passthrough; RT_LOG
      (primitive) non-cancellation: `rde_normal_denominator_field` + `rde_spde_field` +
      `rde_polyrischde_nocancel_field`. Field-correct `k(x)[τ]` algebra (monic-Euclidean
      `rc_gcd`; variable-explicit Quotient/Remainder/ExtendedGCD). Wired into `rt_field_rde`
      log-top ahead of the ansatz; exact tower-identity certified. `Risch`RischDE` extended
      to towers; `Risch`SPDE` + `Risch`PolyRischDENoCancel` boxes exposed. Tests:
      `tests/test_risch_rde_tower.c` (base + SPDE + NoCancel + depth-1/2/3 log towers +
      declines + integrator end-to-end), all green; leak-clean; full transcendental suite +
      P3 suites pass. **NEXT: extract RDE stack → `integrate_risch_rde.{c,h}` (user request).**
- [x] **1b** — RT_EXP non-cancellation: SplitFactor normal/special split in
      `rde_normal_denominator_field`, `RdeSpecialDenomExp` (§6.2 p.190, `ν_τ` valuation
      clears τ-poles), `RdeBoundDegreeExp` (§6.3 p.200); reuses SPDE + NoCancel (deg_τ(b)≥1).
      Closes exp-over-exp Laurent RDEs (`∫ D[E^(E^x)/(1+E^x)] → E^(E^x)/(1+E^x)`), depth-2
      `RischDE[E^(E^x),…]`. All suites green; leak-clean. Deferred to 1c: cancellation
      (`b∈k`, PolyRischDECancelExp) + the ν_τ cancellation refinement (parametric log deriv).
- [ ] **1c** — Cancellation: `PolyRischDECancelPrim` (§6.6 p.212) + `…Exp` (p.213), recurse
      `rde_tower` over `k`. Tests: resonance/coupled-coefficient integrands.
- [ ] **1d** — `rde_weak_normalizer_field` + full `rde_normal_denominator_field` over `k(τ)`
      (spurious-residue strip). Tests: positive-integer-residue tower integrands.
- [ ] **1e** — Delete the `rt_field_rde` `SolveAlways` fallback; general branch becomes
      `rde_tower(i·Dcoef_L, p, ctx@L-1)`. Verify full suite parity + no regressions (A/B any
      pre-existing failures at HEAD).
- [ ] **1f** — Decision wiring: every `rde_tower` NULL → `rt_dec_nonelem`; extend
      `ElementaryIntegralQ` + `strict_unevaluated` guards for the newly-decided class.
- [ ] **(opt) File-split refactor** (Plan §4 Option B): extract RDE stack → `risch_rde.{c,h}`,
      thread `RtDecision*`, add to `tests/CMakeLists.txt`. Behavior-preserving commit.

## Gap 2 — LimitedIntegrate first-class (§7.2)
- [ ] `LimitedIntegrateReduce` (p.248) + §7.1 parametric solve (`LinearConstraints`,
      `ConstantSystem`) in `src/calculus/risch_param.{c,h}`; replace the `SolveAlways`
      approximation in `rt_limited_field_integrate`. Sharpen `RdeBoundDegreePrim` cancellation
      test + `ParametricLogarithmicDerivative` (§7.3). Tests: Bronstein Ex 7.2.x.

## Gap 3 — Tangent tower monomial (`RT_TAN`)
- [ ] `RtKind += RT_TAN`; `rt_tower_build` collects Tan/Cot (special τ²+1), Tanh/Coth (τ²−1).
- [ ] `rde_special_denominator` RT_TAN (`RdeSpecialDenomTan` §6.2 p.192) +
      `rde_polyrischde_cancel` RT_TAN (`PolyRischDECancelTan` §6.6 p.215) → wire
      `Risch\`CoupledDECancelTan`.
- [ ] `rt_field_integrate` RT_TAN dispatch → `IntegrateHypertangent{Reduced,Polynomial}` at
      every level. Tests: nested `Tan[Log[x]]/x`, `1/(a+b Sin[x])`.

## Gap 4 — Retire flat-ansatz cases + remove RT_MAXK
- [ ] After Gaps 1(+3) subsume them: delete/demote `rt_log_tower_case` / `rt_exp_tower_case`,
      remove `nl≤4` caps + `RT_MAXK` (dynamic `RtTower` arrays). Full-suite diff.

## Gap 5 — §5.11 scope boundary (document-only)
- [ ] Note `IntegrateNonLinearNoSpecial` / nonelementary-primitive as out-of-scope-by-
      construction (Bronstein §5.2 p.136) in `RISCH_STATUS.md`; no code.

## Per-increment checklist (CLAUDE.md)
- [ ] `make -j` clean under `-std=c99 -Wall -Wextra`
- [ ] increment tests + full `integrate_risch_transcendental_tests` green
- [ ] valgrind: no new leak blocks over `Sin[1.0]` baseline
- [ ] docs: `docs/spec/builtins/calculus.md` + weekly `docs/spec/changelog/<Mon>.md`;
      `RISCH_STATUS.md` case map; memory `project_risch_p3_decision` linkage
- [ ] commit with Bronstein box→function references

## Review
(fill in as increments land)
