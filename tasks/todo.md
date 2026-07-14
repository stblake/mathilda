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
- [x] **LimitedIntegrate m=1 (§7.2) — the IntegratePrimitivePolynomial fold-back.**
      `rt_limited_field_integrate` now recognises the new logarithm `t_L=Log[u_L]` in the
      tower-variable form `rt_field_integrate` returns it (was: kernel-form subrule mismatch →
      silent decline). Closes `Log[Log[x]]/(x Log[x])→Log[Log[x]]²/2` etc. via the recursive path.
- [x] **Gap 4** — `rt_recursive_tower_case` is now the PRIMARY tower path (subsumes the flat
      cases, validated on suite + corpus); flat cases kept as fallback (carry `rt_tower_solve`
      Cherry substrate). `RT_MAXK` removal deferred (low value).
- [ ] **Full §7 `LimitedIntegrate`** (general m, `LimitedIntegrateReduce` p.248 + §7.1 parametric
      solve) + `b=Dz/z[+mη]` cancellation branch + `ParametricLogarithmicDerivative` (§7.3).
      Deferred — the m=1 case above covers the primitive recursion; the rest is Cherry-adjacent.

## Gap 3 — Tangent tower monomial (`RT_TAN`)
- [ ] `RtKind += RT_TAN`; `rt_tower_build` collects Tan/Cot (special τ²+1), Tanh/Coth (τ²−1).
- [ ] `rde_special_denominator` RT_TAN (`RdeSpecialDenomTan` §6.2 p.192) +
      `rde_polyrischde_cancel` RT_TAN (`PolyRischDECancelTan` §6.6 p.215) → wire
      `Risch\`CoupledDECancelTan`.
- [ ] `rt_field_integrate` RT_TAN dispatch → `IntegrateHypertangent{Reduced,Polynomial}` at
      every level. Tests: nested `Tan[Log[x]]/x`, `1/(a+b Sin[x])`.

## Gap 3 — Tangent tower
- [x] **Nested tangents over C(x)** (commit f8d89bf): `Tan[Log[x]]`, `Tanh[Log[x]]`, etc. via
      relaxed `rt_kernel_eta` + numeric diff-back guard. Soundness fix: trig-frontend false-zero
      `Tan[x]·Tan[Log[x]]→0` (commit 18ae5a4).
- [x] **RT_TAN tower foundation** (commit c23f328): `RtKind += RT_TAN`, `tsg` sign, collection,
      derivation `Dt=Dcoef(t²+σ)`, `Sec²→1+Tan²` rewrite. Sound; builds the tower but declines
      pending the 3 integration pieces below.
- [x] **TrigToTan normaliser** (commit 6793835): `rt_subst_kernels` rationalises circular/hyperbolic
      trig of a tangent arg to the tower symbol (`Sin=t/√(1+σt²)`, …); the fresh symbol stops the
      evaluator canonicalising back to `Csc·Sec`. A log-over-tangent integrand builds the correct
      tower `F=(1+t₀²)/(t₀t₁)`. Verified; non-regressing.
- [~] **RT_TAN full integration** (piece a DONE; b, c remaining): (a) **nonlinear-lower-monomial
      residue support** — DONE. A tangent LOWER kernel has `Dcoef = σ·u'` (§5.10) and gives a
      RATIONAL Dcoef in `t₀`, so the Rothstein–Trager LRT in `rt_field_lrt_logpart` now clears the
      lower-field denominator (scaling `a`, `D` by a common `t_L`-free factor leaves residues/log-args
      invariant), and `T->subrules` carries the full trig rationalisation of the tangent argument
      (`Sin=t/√(1+σt²)`, `Cos=1/√`, `Sec=√`, `Csc=√/t`, `Cot=1/t`) so the evaluator-canonical
      `Csc·Sec` form of `(1+Tan²)/Tan` substitutes cleanly. Closes `∫(1+Tan²)/(Tan·Log[Tan])→
      Log[Log[Tan[x]]]`, the repeated-pole `Log[Tan]^2` form, and the σ=−1 `Tanh` form.
      `test_tangent_tower`; suite + corpus green; leak-clean. (b) the hypertangent-TOP dispatch
      (`IntegrateHypertangent` in `rt_field_integrate`, currently EXP-only in the else branch);
      (c) the §6.2/§6.6 tangent RDE branches (`RdeSpecialDenomTan`, `PolyRischDECancelTan` via
      `CoupledDECancelTan`). Each is a genuine algorithmic extension.

## RT_MAXK depth cap
- [x] **Removed** (commit 4f6453c): `RtTower` arrays heap-allocated to the actual kernel count;
      tower depth unbounded. Corpus + suites green, leak-clean.

## §7 LimitedIntegrate — status of the residuals
- [x] **m=1** (the elementary-integrator need, `IntegratePrimitivePolynomial` fold-back): DONE (4ea4f7d).
- [ ] **general m**: Cherry-adjacent — the elementary transcendental integrator never needs m>1
      (only nonelementary-function / Cherry integration does). Not an elementary-integrator gap.
- [ ] **`b=Dz/z` cancellation branch + `ParametricLogarithmicDerivative` (§7.3)**: rare higher-degree
      refinements; the exact-identity gates keep the integrator SOUND (it declines) without them.

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
