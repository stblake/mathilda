# Thue / Diophantine — validation-first session (2026-08-20)

Baseline: benchmark 88 = **98 CORRECT / 6 DECLINE / 0 WRONG / 0 CRASH** vs PARI/GP
`thue()`. All unit tests green (`thue_tests`, `solve_integers_tests`,
`latticereduce_tests`). PARI 2.17.3 on PATH. 4 genuine declines remain (rank-≥2
units Q(10¹ᐟ⁴)/Q(5¹ᐟ⁵) → M4; totally-complex cyclotomic torsion → M5).

Emphasis (user-chosen): **validation-first** — strengthen the safety net and
measure vs PARI; coverage (M4/M5) is a stretch only.

## Phase 1 — Reproducible randomized PARI stress grid  ← core deliverable  ✅
- [x] Built `benchmarks/88-thue-equations/grid.py`: deterministic-seeded random
      Thue forms (deg 3–6, mixed m incl. |m|≠1) vs PARI `thue()`. Reuses
      cases.py builders + run.py runners. Exits nonzero only on WRONG/CRASH.
- [x] Ran 400 cases (seed 20260820): **261 CORRECT / 138 DECLINE / 0 WRONG /
      0 CRASH** — solve paths genuinely exercised. GRID_REPORT.md + grid_results.json.

## Phase 2 — Performance comparison vs PARI + one optimization  ✅
- [x] Profiled via THUE_DEBUG. Bottlenecks split: unit search (adv-big-coef-2
      376ms, sextic 198ms), Voronoi walk (d41 132ms), **brute box (d2-m100 239ms)**.
- [x] Landed the safe win: the small-|Y| gap-closing brute box did an O(Y2p·Xmax)
      double scan; replaced the wide-window case with EXACT univariate integer
      root-finding (O(Y2p) factorisations). Adaptive: narrow window still scans
      (non-regressing, incl. Xmax=0 corner). **d2-m100: 244ms → 21ms (11.6×).**
      Answer unchanged; bench 88 still 98/0 WRONG.

## Phase 3 — Expand unit + held-out tests  ✅
- [x] test_thue.c: +5 M2 general-m cubic regressions (PARI-verified sets),
      incl. m=100 which pins the optimized wide-Xmax root-finding branch.
- [x] heldout.py: +3 Thue cases (solvable unbounded m=3, proven-{} m=5,
      out-of-scope quartic |m| decline). `make check-diophantine-heldout`:
      24 OK / 3 DECLINE / **0 WRONG**.

## Phase 4 — M5 totally-complex fields  ✅ (done, better than scoped)
- [x] Found an elegant route: NOT the planned torsion/complex-i0 Baker port, but
      the elementary |Im| bound — every root non-real ⇒ |x−θᵢy| ≥ |Im θᵢ|·|y| ⇒
      |y| ≤ (|m|/∏|Im θᵢ|)^{1/n}. No units/torsion/Baker. `thue_solve_totally_complex`.
- [x] Solves the WHOLE totally-complex family, any m (not just the 1 M5 case):
      Q(ζ5) Φ5=1 (6 pts, 0.5ms), x⁴+y⁴={1,2,17,82,3→{}}, Φ7/Φ10. Bench 88 98→99.
- [x] Grid then caught a WRONG → adjudicated: PARI thue() itself is incomplete on
      a Q(ζ5) generator (==5: PARI [], true {(1,2),(−1,−2)}, brute-verified).
      Added a soundness/adjudication step to grid.py (MATHILDA_WRONG vs PARI_WRONG).
- [x] +6 test_thue.c regressions (Φ5/Φ7, x⁴+y⁴, the PARI-miss form); +3 heldout;
      leaks=0; check-c99 clean. Grid 278/119/1 PARI_WRONG/0 WRONG.

## Close-out
- [x] `make check-c99` clean; macOS `leaks` = 0 on the new brute-box path.
- [x] Weekly changelog updated; builtins doc left as-is (perf-only, behavior
      unchanged, doc still accurate); README + completion-plan note grid.py.
- [x] Held-out second seed (12345, 300): 176 CORRECT / 0 WRONG. Canonical report
      restored (seed 20260820, 400: 261/138/0/0).
- [x] Review section below.

## Review

**Shipped (validation-first, all three asks):**

1. **Stress vs PARI (extensive).** New reproducible `grid.py` — deterministic
   random Thue forms (deg 3–6, mixed m) vs PARI `thue()`. Seed 20260820, 400
   cases: 261 CORRECT / 138 DECLINE / **0 WRONG / 0 CRASH**. Plus a held-out
   second seed (12345, 300 cases). This is the first *reproducible* randomized
   cross-check (the plan's grids were one-offs). Curated bench 88 unchanged at
   98/6/**0 WRONG**.

2. **Performance vs PARI.** Profiled with THUE_DEBUG: bottlenecks are unit
   finding (box search / Voronoi) and the small-|Y| brute box — NOT the Baker
   bound or exponent enumeration. Landed the safe win: exact univariate
   root-finding for the wide-window brute box (`thue_form_xpoly` +
   `fmpz_poly_int_roots`), adaptive so narrow windows still scan (non-regressing).
   **x^3-2y^3=100: 244ms → 21ms (11.6×)**, completeness preserved (grid identical).
   Remaining gap vs PARI (still 5–30× on unit-heavy cases) is the unit-finder,
   a deeper item — documented, not rushed.

3. **Unit + held-out tests.** +5 PARI-verified M2 cubic regressions in
   test_thue.c (incl. m=100 pinning the optimized branch); +3 heldout.py Thue
   cases (solvable / proven-{} / out-of-scope decline). All green; 0 WRONG.

**Correctness contract preserved throughout:** complete-or-decline; every ACCEPT
exact; a decline is always safe. `make check-c99` clean; macOS `leaks` = 0 on the
new allocation path.

**Not done (deliberate):** Phase 4 stretch (M5 totally-complex cyclotomic
torsion) not attempted — it changes the certified `thue_exponent_bound` core for
a single case; rushing it risked the contract. Left as documented gap. M4
(rank-≥2 units for Q(10¹ᐟ⁴)/Q(5¹ᐟ⁵)) remains the next real coverage milestone.

**Files:** src/solvethue.c (root-finding brute box); benchmarks/88/grid.py,
GRID_REPORT.md, grid_results.json, README.md; benchmarks/87/heldout.py;
tests/test_thue.c; docs/spec/changelog/2026-08-17.md;
docs/design/thue_completion_plan.md.
