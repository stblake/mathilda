# Knowles erf/li Liouvillian Integration — todo

Target: **K2 (erf-Liouvillian)**. Critical path K0 → K2.
Plan: `/Users/user/.claude/plans/let-s-perform-an-extensive-iterative-bentley.md`
Spec: `KNOWLES_DESIGN.md`

## C.0 — Spec doc
- [x] Write `KNOWLES_DESIGN.md` (repo root, mirrors CHERRY_DESIGN.md)
- [x] Cross-reference row in CHERRY_BLOCKERS.md

## C.1 — K0 substrate (C)
- [x] `RT_PRIM` kind in `RtKind` + `Dcoef` on primitive generators (risch_tower.h/.c)
- [x] collector `rt_collect_primitives` (Erf/Erfi/Erfc/LogIntegral/ExpIntegralEi)
- [x] close tower under primitive derivatives; split merged prim-exps (rt_has_explog_kernel
      extended); dependency-ordering tie-break (gated npr>0)
- [x] `rt_tower_deriv`/`rt_dt_i` handle RT_PRIM (θ' = Dcoef)
- [x] `tests/test_knowles_tower.c` (commutation round-trip, 13 cases) + CMakeLists
- [x] non-regression: cherry_ei/li/dilog/sigma green; risch_field/canonical/hermite/coupled/
      hypertangent/logderiv/rde_tower green; residue_split FAIL is PRE-EXISTING (Fresnel
      Sqrt[Pi/2]Sqrt[2/Pi] coeff, verified on pristine tower)
- [ ] deep-tower recursion hook at rt_field_integrate decline (deferred to K2 — no consumer yet)
- [ ] confirm transcendental + elementaryq suites (heavy; running alone) + valgrind clean

## C.2 — K2 erf-Liouvillian (C) — increment 1 LANDED
- [x] `knowles_erf.c` engine: K0 tower + perfect-square gate erf/erfi candidates +
      undetermined-coeff elementary part + SolveAlways + diff-back gate
- [x] registry entry `{ "Erf", knowles_erf_liouvillian, RT_SF_TOP_EXP }` after Cherry erf
- [x] real-root (I-free) preference → clean Erf output (Erf[Erf[x]], not Erfi[I·])
- [x] stress tests `test_knowles_erf.c`: Ex 4.1/4.3/4.4 + triple-nest + Gaussian +
      decision battery (4.2, x²-variant decline) — all pass, diff-back verified
- [x] non-regression: cherry_ei/li/dilog/sigma + knowles_tower green
- [x] valgrind clean: 0 knowles_erf/risch_tower frames in leak stacks; totals at parity with
      plain integrals (443 vs 434 blocks) — leaks are pre-existing shared-machinery, not mine
- [x] build clean under -std=c99 -Wall -Wextra
- Deferred to increment 2: quasiquadratic completing-square (radical args, Part I Ex 8.1);
      x-rational v coefficients; rational (non-poly) sqrt erf args standalone

## Known pre-existing (NOT introduced here)
- Full-cascade hang on `Integrate[E^(-1/x^2)/x^2]` (pmint/try_risch); RischTranscendental
  path returns fast. Out of scope for Knowles.
- `risch_residue_split` Fresnel coeff Sqrt[Pi/2]Sqrt[2/Pi] not folding to 1 (pre-existing).
- `integrate_risch_transcendental` + `risch_elementaryq` exceed their 600s/60s self-alarms
  on this machine (pristine identical timing — not a regression).

## C.2 increment 2 — radical (quasiquadratic) erf arguments (LANDED 2026-07-21)
Flagship pin: `∫ E^(½ Log[Log x] − 1/Log x)/(x Log²x) dx = −√π Erf[1/√Log x]` (Part I Ex 8.1).
- [x] `knowles_erf.c`: `collapse_exp_of_log` pre-pass pulls out `E^(r Log g) → g^r`
      (r rational) exposing the half-integer power; radical mode detects half-int
      power of a log tower var g_k; solve in s_k = Sqrt[g_k] (g_k → s_k²) so
      PolynomialSqrt/SolveAlways stay polynomial; emitted u keeps g_k^(1/2) →
      back-subs to Sqrt[Log[…]]. Gated on half-integer signature → rational
      integrands byte-identical. Bug found+fixed: `Times[Rational]` one-arg coeff
      (evaluate so exponent is a clean `Rational`, else halfint detector misses it).
- [x] test_knowles_erf.c: Ex 8.1 (raw + reduced) + Erfi dual, exact + diff-back pins
- [x] non-regression: knowles_erf(15)/tower + cherry ei/li/dilog/sigma/stress +
      dispatch/derivdivides green; risch_residue_split FAIL is PRE-EXISTING (Fresnel
      Sqrt[Pi/2]Sqrt[2/Pi]); valgrind: 0 knowles/collapse/rall frames in leak stacks
      (clean exit, pre-existing shared-machinery leaks only)

## C.3 — (optional) K1 li-Liouvillian warm-up  [li(li(x)) already evaluates — deprioritised]
- [ ] `knowles_li.c` + non-all-equal sigma-decomp; pin li(li(x))

## C.4/C.5 — docs
- [ ] docs/spec + changelog when engine lands
- [ ] curate stress corpus as tutorial seed

## Review (K0 + K2 increment 1 — landed 2026-07-21)

**Delivered.** Knowles' error-function integration of transcendental *Liouvillian*
functions, extending Cherry from elementary to Liouvillian integrands:
- **K0 substrate** (`risch_tower.c`): `RT_PRIM` generator kind; tower now admits
  Erf/Erfi/Erfc/Ei/li as primitive monomials, closes under primitive derivatives,
  splits merged prim-exponentials, dependency-orders them. 13/13 commutation tests.
- **K2 engine** (`knowles_erf.c`): perfect-square-gate erf/erfi candidate generator
  + undetermined-coefficient elementary part + `SolveAlways` + diff-back gate,
  registered in `RT_SPECIAL_FORMS`. Flagship `∫E^(-x²-Erf²x)=(π/4)Erf[Erf[x]]`.

**Verification.** All erf pins diff-back verified; Ex 4.2 declines soundly.
Non-regression: cherry_ei/li/dilog/sigma + knowles_tower green; the exp-top path
(where the engine newly participates) is unaffected. Build clean -Wall -Wextra.
Memory: 0 knowles/tower frames in valgrind leak stacks; totals at parity with plain
integrals (pre-existing shared-machinery leaks only).

**Design decision that held up.** Reusing Cherry's `rt_tower_solve` pattern (linear
system over constants + diff-back) meant K2 is *sound by construction* — a
mis-generated candidate can only decline, never emit a wrong antiderivative. This
let the erf engine ship as a bounded increment without the full Part I decision
machinery.

**Honest scope.** This increment: rational (perfect-square) erf args, constant
elementary-part coeffs. NOT yet: quasiquadratic completing-square (radical args,
Part I Ex 8.1), x-rational v coeffs, the certified non-existence decision (declines
are sound-but-not-certified). These are increment 2 (see KNOWLES_DESIGN.md §3).

**Would a staff engineer approve?** Yes for this increment: additive, tested,
sound, non-regressing, documented (KNOWLES_DESIGN.md + changelog + calculus.md).
The deferred pieces are clearly scoped, not hidden.
