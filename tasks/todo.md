# M45 — §2.2.26 corpus (Problems 2501–2600), push toward the honest ceiling

## Plan
Add Nasser Abbasi's §2.2.26 (2501–2600) to the DSolve corpus; measure; close
honest gaps with **general root-cause fixes only**. Section: 100 scalar ODEs
(50 IVP, 0 systems) — first-order-nonlinear-heavy (homogeneous class A/C,
dAlembert, Bernoulli, Abel, exact, separable, Riccati) + a 2nd-order block
(const-coeff "missing x", Euler–Cauchy, Emden–Fowler, Gegenbauer, with-symmetry).
Invariants: **0 FAIL (0 wrong answers), 0 regression**; real methods only — a
research-grade case is a *documented decline* (all residue `sympySolved=False`),
never a faked pass. Leave uncommitted. Milestone M45, version 0.142 → 0.143.
Planning baseline (Sep-13 binary): **87/100, 0 FAIL, 0 crash, 13 UNEVAL**.
Targets: **2532** (Bernoulli `DSolve\`Y` leak, sympy=True), **2592** (2nd-order
linear rational-coeff — same class as §2.2.25-2410). Target ~89/100.

## Tasks
### Stage 0 — ingest + measure
- [x] 1. Generated `DSolve_test_status/DE_examples_2226.m` (100 recs, 50 IVP, 0 sys)
- [x] 2. Converter output vetted (t/x split legit, params OK, no glued symbols, no
      t/x-mixed eq). No converter fix needed → §2.2.20–25 byte-identical (converter untouched).
- [x] 3. Fresh rebuild + `dsolve_corpus_tests`; baseline **87/100, 0 FAIL, 0 crash, 13 UNEVAL**.
- [x] 4. Generated `reports/2.2.26.{md,tsv}` (+ regenerated `2.2.25.{md,tsv}`).
- [x] 5. Triaged: 2532 (Bernoulli gap, sympy=True), 2592 (poly-sol gap), 11 honest declines.

### Stage 1 — root-cause gap closure (0 FAIL, 0 regression)
- [x] 6. 0 FAIL confirmed on fresh binary.
- [x] 7. **Fix A — Bernoulli `DSolve\`Y` leak** (`dsolve_bernoulli.c`): `Cancel` the linearised
      coefficients before the linear solve; 2532 → clean integral-form → **PASS**.
- [x] 8. **Fix B — polynomial-solution 2nd-order** (new `dsolve_ratsol2.c`, before Kovacic;
      self-contained VoP via `dsolve_variation_of_parameters`): **2592 + §2.2.25-2410 → PASS**.
- [x] 9. Residue: §2.2.26 = 11, §2.2.25 = 3 — ALL `sympySolved=False` (honest).

### Stage 2 — register + docs + version
- [x] 10. Saved `reports/2.2.26.{tsv,md}` + regen `2.2.25.{tsv,md}`; gate `dsolve_corpus_2_2_26_tests`
      (baseline 11); §2.2.25 gate 4→3.
- [x] 11. STATUS.md §2.2.26 block + §2.2.25 M45 row; README.md rows.
- [x] 12. `src/version.h` 0.142→0.143; DSOLVE_PLAN.md M45 entry; new changelog
      `docs/spec/changelog/2026-09-14.md` + Mathilda_spec.md table row.

### Stage 3 — verification
- [x] 13. `ctest -R dsolve_corpus_2_2`: §2.2.26 (89, baseline 11), §2.2.25 (97, baseline 3),
      and all other sections PASS **in isolation**. Full back-to-back sweep flaked 3 sections
      (2.2.12/15/20) — PROVEN not an M45 regression: every flaky case runs byte-identical code
      (1st-order, systems, high-order → ratsol2 guard-declines at n!=2/nfun!=1; Bernoulli Cancel
      only on the Bernoulli success path). Pre-existing 8s-`TimeConstrained` boundary flakiness
      under sweep CPU contention (local-only gates; CI runs only check-c99 + Linux build).
- [x] 14a. `make check-c99` PASSES. REPL spot-checks: 2532/2592/2410 clean, no `DSolve\`Y`;
      regressions (const-coeff, Euler, y''=0) unchanged.
- [x] 14b. Optimized ratsol2 decline path (cap 10→6, CoefficientList not Coefficient-loop):
      §2.2.20-1960 (Frobenius, ratsol2 declines) 6.8s→5.7s — MORE margin, not less.
- [x] 15. valgrind leak-neutral: definite-loss (~13.5–14 KB, equation-dependent) is the
      pre-existing per-call Integrate/Simplify leak reached via shared helpers
      (simp_canon.c, dsolve_linear_factor_solve); ZERO leaked allocation in my ratsol2 frames
      or the Bernoulli Cancel lines. `dsolve_tests` SIGALRM 142 (Risch-Norman enum cap) is
      pre-existing and untouched by M45 (no Risch-Norman code changed).
- [x] 16. Staff-engineer self-review of `dsolve_ratsol2.c` memory: clean on all 5 exit paths.

## Review

**Outcome.** M45 adds §2.2.26 (Problems 2501–2600) to the DSolve corpus at **89/100, 0 FAIL,
0 crash**, and lifts §2.2.25 **96 → 97/100** as a bonus, with **0 regression** (every §2.2.x
section passes its gate in isolation). Version 0.142 → 0.143.

**Root-cause fixes (general, no overfit):**
- **Bernoulli `DSolve\`Y` integrating-factor leak** (`dsolve_bernoulli.c`): `Cancel` the
  linearised coefficients to lowest terms in the reduction variable before the linear solve,
  so the frozen internal symbol never leaks into a non-elementary `Integrate`. Fixes 2532 and
  the whole non-elementary-IF Bernoulli class; Y-free coefficients unchanged.
- **New `DSolve\`PolynomialSolution` method** (`dsolve_ratsol2.c`, before Kovacic): degree-bounded
  polynomial-fundamental-set search + VoP over the clean basis. Fixes 2592 AND §2.2.25-2410 —
  the class Kovacic's complex-radical basis (√(t−i)√(t+i)) blocked. Self-verifying; cheap decline.

**Honest ceiling.** §2.2.26 residue 11 and §2.2.25 residue 3 are ALL `sympySolved=False`
(non-integrable Riccati/Abel, implicit y=G/x=G forms, abstract-coefficient, parabolic-cylinder,
non-elementary IF) — genuinely beyond current CAS reach, documented declines, no wrong answers.

**Deliverables.** New `DE_examples_2226.m` + `reports/2.2.26.{md,tsv}`; regenerated
`reports/2.2.25.{md,tsv}`; gate `dsolve_corpus_2_2_26_tests` (baseline 11); §2.2.25 gate 4→3;
new `src/calculus/dsolve_ratsol2.c` (+ dsolve.c wiring, tests/CMakeLists.txt);
STATUS.md / README.md / DSOLVE_PLAN.md / new weekly changelog `2026-09-14.md` / Mathilda_spec.md;
version.h 0.143. Left uncommitted.

**Pre-existing (not addressed, not caused here):** the ~13.5 KB per-call Integrate/Simplify
leak; `dsolve_tests` SIGALRM 142 (Risch-Norman); full-sweep 8s-boundary flakiness on a few
1st-order/system/high-order cases in §2.2.12/15/20 (each passes in isolation).
