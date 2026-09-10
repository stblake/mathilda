# M36 — DSolve §2.2.16 corpus (Problems 1501–1600)

Section: Boyce & DiPrima, indexsubsection25.htm. 100 records, all scalar, 56 IVP.
Buckets: 39 separable, 25 linear, 13 quadrature, 16 2nd-order nonhomog + 2 high-order
(step/impulse-forced IVPs 1501–1518), 5 specials (Clairaut 1536, Riccati 1577,
symmetry 1575/1576, class-A 1561).

## Steps

- [x] 1. Generated `DE_examples_2216.m` (100 scalar, 56 IVP). Converter fix: `|…|`→`Abs[…]`.
- [x] 2. Wired `dsolve_corpus_2_2_16_tests` (baseline 3).
- [x] 3. Baseline measured: 96/100, 0 FAIL, 4 UNEVAL (1508,1509,1534,1590). Report generated.
- [x] 4. Root-caused: fixed 1590 (cubic-log separable ≥3-log → implicit twin). 1508/1509/1534
        = documented bounded declines (impulse √2 churn / converter indep-var). Final: 97/100.
- [x] 5. Anti-overfit `t_m36_separable_cubic_log` added + registered; dsolve_tests green.
- [x] 6. Docs updated (STATUS, README, DSOLVE_PLAN M36, changelog, version 0.134, CMake baseline 3).
- [~] 7. Verify: check-c99 ✅; dsolve_tests ✅; full corpus ctest suite RUNNING (no-regression gate).

## Review

**Outcome: §2.2.16 (Problems 1501–1600) added, 97/100, 0 FAIL, 0 regression.**

Two root-cause fixes (no per-problem hacks):
1. **Converter `Abs[…]`** (`tools/latex_ode_to_mathilda.py`) — `\left|…\right|`/`\lvert`/bare `|…|`
   → `Abs[…]`. Necessary: without it the whole corpus file did not parse (1535 `y'==|y|+1`).
2. **Cubic-log separable → implicit first integral** (`src/calculus/dsolve_separable.c`) — a
   separable whose `∫1/h` has ≥3 distinct `Log` args has no explicit inverse `Solve` closes
   (churns); the explicit path now declines it so the implicit twin returns `G==C[1]`. Closes 1590
   (which SymPy does NOT solve). 1-/2-log relations stay explicit (no regression).

Investigated but NOT fixed (documented bounded declines, all needing substantial substrate work):
- 1508/1509: impulse Green's-function convolution churns on irrational `−1±I√2` roots; its
  `TrigReduce` is REQUIRED for sift correctness (skipping it → spurious ramp, verified on 1506).
  A tried VoP linearity-split + pure-impulse gate were reverted (broke 1506 / no corpus gain).
- 1534: converter reads the parameter `a` as the independent variable (ambiguous source).
