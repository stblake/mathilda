# Task: Simplify over Root objects — stop the hang + fully-general algebraic collapse

Plan: /Users/user/.claude/plans/the-following-hangs-in-declarative-cascade.md
Scope: fully general reach-0.

## Milestone 0 — reproduce & baseline  ✅
- [x] Build `./Mathilda` (v0.164, current)
- [x] Reproduced hang: `Simplify[d]` EXIT=124 at 45s; PZQ True; LeafCount 192
- [x] Bare-Root: Simplify leaves `Root^3+Root` inert; `RootReduce[Root^3+Root]`→I√2

## Milestone 1 — never hangs (Part B)  ✅
- [x] Promote `simp_contains_root_head` to shared helper (simp_internal.h)
- [x] Gate polish Factor (simp_builtins.c:1083) with `!simp_contains_root_head`
- [x] Gate alg-Plus Together (simp_builtins.c:899)
- [x] **ROOT-CAUSE FIX**: `collect_variables` (poly.c:575) treats Root as a
      constant, not a variable. THIS is what killed the hang (the two gates alone
      did NOT — the entry was the round-loop/bottom-up together/cancel recursion).
- [x] Harden `exact_poly_div` (poly.c:1415) with `prev_degR` guard (defense-in-depth)
- [x] Verify: `Simplify[d]` now returns (EXIT=0); result terminates, not yet 0
- [x] Sanity: Variables[Root+x]={x}; non-Root Factor/Together/Cancel/GCD/Solve OK
- [ ] Regression suite (building in bg)
- Note: transform_can_fire Together/Cancel backstop NOT needed (collect_variables
  fix is the real cure). Skipped to avoid over-restriction.

## Milestone 2 — Root arithmetic (Parts A + D)  ✅
- [x] Promote `rr_thread_coeffs` → `flint_qqbar_reduce_coeffs` (flint_qqbar.c/.h)
- [x] Route rootreduce.c through it; delete static copy; RootReduce tests pass
- [x] Add Phase-0d coefficient pre-pass (simp_search.c) gated on Root head.
      (Skipped relaxing the RootReduce gate — Phase 0d handles bare Root arith.)
- [x] Verify: `Root^3 + Root` → I√2; `Root-Root` → 0

## Milestone 3 — shown example → 0 (Part C, trigrat)  ✅
- [x] NOTE: collect-by-generators + qqbar-reduce theory was WRONG (numerator is
      zero only mod the trig ideal, not as a free polynomial). Discarded.
- [x] CORRECT FIX: trigrat wrapper converts constant-algebraic Roots to radicals
      (ToRadicals) at entry; existing radical machinery (Sqrt[2]-proven) reaches 0.
      Degree>=5 non-radical Roots left untouched (graceful, no hang).
- [x] TRIGRAT_OPAQUE_MAX=24 already >= 4 (no change needed)
- [x] Verify: `Simplify[d]` → 0 (8.6ms); `D[F,x]//Simplify` → Sqrt[Tan[x]] (9.6ms)

## Milestone 4 — fully general (non-trig)  ✅
- [x] Works via collect_variables(Root=constant) + Together-groups-by-x + Phase-0d
      qqbar coeff reduction. No separate front-end needed.
- [x] Verify: `D[Integrate[x/(x^3-2),x]] - x/(x^3-2)` → 0 (66ms)

## Closeout (Part E)  ✅
- [x] Tests: test_simplify_hang.c (forked termination case + 6 exact assertions)
- [x] `make check-c99` passes; macOS `leaks` on new paths = 0 leaks
- [x] Bump src/version.h → 0.165; docs/spec/builtins/simplification.md + changelog
- [x] memory: project_simplify_root_objects.md (+ MEMORY.md index)
- [ ] Commit + tag v0.165 — LEFT TO USER (commit only when asked)

## Review

**Delivered (v0.165):** Simplify over `Root[...]` objects — robust + fully-general
reach-0.
1. **Root cause / never hangs:** `collect_variables` (poly.c) treats a Root object
   as a constant, not a polynomial variable — kills the multivariate blow-up at
   source for every poly-engine caller. + polish-Factor/alg-Together gates +
   `exact_poly_div` degree guard.
2. **Root arithmetic:** `flint_qqbar_reduce_coeffs` (promoted from rr_thread_coeffs)
   + Phase-0d Simplify pre-pass → `Root^3+Root` → `I√2`.
3. **Reach 0:** trigrat ToRadicals-at-entry for radical-expressible Roots →
   `Simplify[D[∫√(tan x)]−√(tan x)]` → 0 (12 ms, was >300 s hang); `D[F]//Simplify`
   → `Sqrt[Tan[x]]`. Non-trig identities collapse via Together + qqbar pass.

**Wrong turn corrected:** first tried collect-by-{s,c,l} + qqbar-reduce coefficients
— WRONG (numerator is zero only mod the trig ideal). ToRadicals→existing machinery is
correct and simpler.

**Validation:** 50+ targeted suites green (poly/rat/factor/mvfactor/simplify/
simplify_hang/trigrat/rootreduce/radical/solve/algebraicnumber/minpoly/qafactor/
fullsimplify/integrate/…). Pre-existing failures (risch_rde_tower, dsolve_stress,
dsolve_corpus, moebiusmu/primenu, qrdecomposition segfaults) confirmed identical on a
stash-reverted clean baseline — NOT caused by this change.

**Limit (honest):** a genuinely non-radical (degree ≥5) Root coefficient does not
collapse to 0, but degrades gracefully (no hang, unchanged form). Covering it would
need common-field AlgebraicNumber arithmetic — deferred.
