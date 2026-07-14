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
- [x] **1c** — Cancellation (`b∈k*`, deg_τ(b)=0): `rde_polyrischde_cancel`
      (`PolyRischDECancel{Prim,Exp}`, §6.6 p.212/213) builds q top-down, per-coefficient
      lower-field Risch DE `D[s]+(b+m·η)s=lc_τ(c)` via **recursive `rde_tower` over `K_{m-1}`**.
      Closes const-coeff RDEs (`[1,1/x+Log[x]]→Log[x]`, `[2,E^(-x)]→E^(-x)`), cancel→nocancel
      chains (`1/(1+E^(E^x))`), integrator `∫D[E^x/(1+Log[x])]`. All suites green; leak-clean.
      Deferred to Gap 2: the `b=Dz/z[+mη]` antidifferentiation branch (LimitedIntegrate).
- [~] **1d** — `rde_weak_normalizer_field`: **DEFERRED — no reachable test case.** Every
      constructible positive-integer-residue tower RDE already solves without it
      (`RdeNormalDenominator` h + SPDE + cancel find them; exact-identity gate keeps it
      sound). Theoretically-complete-but-pre-empted, like the ansatz resonance code. Revisit
      only if a real integrand is found that needs it.
- [x] **1e** — Retired the `rt_field_rde` `SolveAlways` `h/pd` ansatz (~158 lines) + orphaned
      `rt_resonance_int`; the general branch routes every field RDE through `rde_tower` and its
      NULL is an authoritative "no rational solution in K_L". Also extended the rde_tower gate to
      RT_EXP tops. Verified non-regressing: 294-assertion transcendental suite + broad
      `integrals_tests` corpus + all ansatz-era field-RDE examples still close. Leak-clean.
- [x] **1f** — Decision wiring inherent in 1e: `rt_field_rde` NULL → `rt_dec_nonelem` (matches
      the ansatz's prior authoritative-NULL behavior); `ElementaryIntegralQ` suite green. Residual
      non-authoritative declines (tangent top, `b=Dz/z`) are out of tower scope / rare.
- [ ] **(opt) File-split refactor** (Plan §4 Option B): extract RDE stack → `risch_rde.{c,h}`,
      thread `RtDecision*`, add to `tests/CMakeLists.txt`. Behavior-preserving commit.

## Gap 2 — antidifferentiation / LimitedIntegrate
- [x] **Antidifferentiation branch** (`sp.b=0`, `D h = c`): `rde_tower` integrates `c` in
      `K_m` via `rt_field_integrate` (RDE↔integrator mutual recursion), gated on a rational
      tower result (`rc_is_tower_rational`; a new log ⇒ decline). Hardened the identity gate
      (`Expand` the residual numerator — `Together` doesn't expand products). Closes
      `RischDE[0, x E^x]→(x−1)E^x`, `RischDE[0, Log[x]]→x Log[x]−x`, `∫E^(E^x)(x E^x+1)→x E^(E^x)`.
      Tests: `test_tower_antidiff`. Leak-clean; all suites green.
- [ ] **Full §7 `LimitedIntegrate`** (`LimitedIntegrateReduce` p.248 + §7.1 parametric solve)
      as a first-class facility in `src/calculus/risch_param.{c,h}`, replacing the `SolveAlways`
      approximation in `rt_limited_field_integrate`; the `b = Dz/z [+ mη]` cancellation branch;
      `ParametricLogarithmicDerivative` (§7.3). Deferred — the integrator already approximates.

## Gap 3 — Tangent tower monomial (`RT_TAN`)
- [ ] `RtKind += RT_TAN`; `rt_tower_build` collects Tan/Cot (special τ²+1), Tanh/Coth (τ²−1).
- [ ] `rde_special_denominator` RT_TAN (`RdeSpecialDenomTan` §6.2 p.192) +
      `rde_polyrischde_cancel` RT_TAN (`PolyRischDECancelTan` §6.6 p.215) → wire
      `Risch\`CoupledDECancelTan`.
- [ ] `rt_field_integrate` RT_TAN dispatch → `IntegrateHypertangent{Reduced,Polynomial}` at
      every level. Tests: nested `Tan[Log[x]]/x`, `1/(a+b Sin[x])`.

## Gap 3 — Tangent tower
- [x] **Nested tangents over C(x)** (commit f8d89bf): `Tan[Log[x]]`, `Tanh[Log[x]]`, etc. via
      relaxed `rt_kernel_eta` (eta-kernel-free = genuine over-C(x)) + numeric diff-back guard.
      Plus a soundness fix: trig-frontend false-zero `Tan[x]·Tan[Log[x]]→0` (commit 18ae5a4).
- [ ] **Full `RT_TAN` tower monomial** — tangent MIXED with an independent Log/Exp of x
      (`Tan[x]·Log[x]`, `Log[Tan[x]]`) + the §6.2/§6.6 tangent RDE branches (`RdeSpecialDenomTan`,
      `PolyRischDECancelTan`) wiring `CoupledDECancelTan`. Large structural build; deferred.

## Gap 4 — Retire flat-ansatz cases + remove RT_MAXK
- [~] **BLOCKED (experiment done).** Disabling `rt_log_tower_case`/`rt_exp_tower_case` regresses
      the primitive-polynomial class `Log[Log[x]]/(x Log[x]) → Log[Log[x]]²/2`,
      `Log[Log[x]]^5/(x Log[x])`: the recursive path can't fold back the new logarithm (needs
      §5.8 `IntegratePrimitivePolynomial` / `LimitedIntegrate`, the same §7 residual). The flat
      cases' SolveAlways ansatz finds it. **Prerequisite: full `LimitedIntegrate`.**
- [ ] `RT_MAXK`=5 depth-cap removal (dynamic `RtTower` arrays) — low value (depth-5 already deep).

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
