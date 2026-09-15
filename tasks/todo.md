# DSolve M47 — §2.2.28 corpus (Problems 2701–2800)

Milestone: M47 = §2.2.28. Version bump 0.144 → 0.145.
Definition of done: 0 FAIL, 0 crash; UNEVAL OK only for no-closed-form residue.

## Phase 1 — Fetch & convert
- [x] Fetch Ch2.S2.SS28.htm (browser UA), verified title "Problems 2701 to 2800"
- [x] Convert → DE_examples_2228.m (100 records: 18 scalar [3 IVP], 82 systems)
- [x] Converter fix: subscripted arbitrary funcs f_1(t)/f_2(t) → f1[t]/f2[t]
      (2708/2762/2780); byte-identical regen of §2.2.27 confirmed

## Phase 2 — Baseline measurement
- [x] Build dsolve_corpus_tests
- [x] Baseline: 84/100 PASS, 0 FAIL, 0 crash, 16 UNEVAL
- [x] reports/2.2.28.md + .tsv generated

## Phase 3 — Root-cause fixes (2 general fixes, 0 FAIL throughout)
- [x] 2719 (y''''+y==g[t]) radical-root VoP integrand simplify
      (ds_has_var_fractional_power guard; 2406 t^(5/2) hang untouched)
- [x] 2713 zero-ladder homogeneous IVP → y≡0 (tight guard; BVP excluded)
- [x] Re-run §2.2.28: 86/100, 0 FAIL, scalars 18/18
- [x] Re-run all gates → 0 regression (100% pass, 28/28; §2.1.2 no-IC + VoP-validated)

## Phase 4 — Record & land
- [x] tests/CMakeLists.txt — dsolve_corpus_2_2_28_tests gate (baseline 14)
- [x] STATUS.md — §2.2.28 block + M47 wave-history line
- [x] README.md — DE_examples_2228.m row
- [x] DSOLVE_PLAN.md — M47 bullet
- [x] docs/spec/changelog/2026-09-14.md — M47 entry
- [x] src/version.h — 0.144 → 0.145; rebuilt (banner shows 0.145)
- [x] make check-c99 clean
- [x] Confirm regression run: all gates hold (100% pass, 0 FAIL)
- [x] Confirm §2.2.28 ctest gate passes at baseline 14 (Passed 70.8s)
- [x] Review section (below)

## Residue 14 (all systems, no scalar gap)
- 8 nonlinear systems 2788–2795 (sympy=False, no closed form)
- 3×3 E^t system 2785 (sympy=False)
- system arbitrary forcing 2708/2762/2780 (needs system-level VoP)
- system DiracDelta/UnitStep 2781/2782 (needs systems Laplace/Green's)

## Review

**Result: §2.2.28 landed at 86/100, 0 FAIL, 0 crash (scalars 18/18). v0.145.**

Two general root-cause fixes (both in `src/calculus/dsolve_common.c`):
1. Constant-base-radical VoP integrand simplify (`dsolve_variation_of_parameters`) —
   fixed 2719 (`y''''+y==g[t]`, irrational roots `(±1±i)/√2`); new
   `ds_has_var_fractional_power` keeps the §2.2.25-2406 `t^(5/2)` Simplify hang guarded off.
2. Homogeneous-linear-IVP zero-ladder → `y≡0` (`dsolve_fit_constants`) — fixed 2713
   (irreducible-quartic `Root[]` fundamental set). Guarded by `ds_is_structural_zero(body[C->0])`
   so nonhomogeneous IVPs fall through.

Converter fix (`tools/latex_ode_to_mathilda.py`): subscripted arbitrary forcing functions
`f_1(t)`/`f_2(t)` → `f1[t]`/`f2[t]` (2708/2762/2780); §2.2.27 regenerates byte-identically.

**Process note (caught + fixed in-session):** the first cut of the zero-IVP shortcut lacked
the homogeneity check and shipped WRONG answers (15 sections FAILed in the full regression —
`y''+y==t, zero ICs` → particular `t`, which fails `y'(0)=0`). The `ctest -R dsolve_corpus`
regression caught it; fixed with `ds_is_structural_zero(body[C->0])`, re-run 100% green.
Lesson recorded in `tasks/lessons.md` + memory. Reinforces: always run the FULL corpus
regression before landing a change to a shared path, even when the target section is green.

**Verification:** `make check-c99` clean; `./Mathilda -v` → 0.145; `ctest -R dsolve_corpus`
100% (28/28) incl. §2.1.2. Not committed (awaiting user go-ahead; would branch off main first).
