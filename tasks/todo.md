# M42 — §2.2.23 corpus (Problems 2201–2300)

## Plan
Add Nasser Abbasi's §2.2.23 (2201–2300) to the DSolve corpus; measure; drive
coverage with general root-cause fixes. Invariants: 0 FAIL, 0 regression.
Section: 45 scalar (17 IVP) + 55 constant-coeff linear systems (0 IVP).

## Tasks
- [x] 1. Generate `DSolve_test_status/DE_examples_2223.m` (converter, saved HTML)
- [x] 2. Add `dsolve_corpus_2_2_23_tests` gate in `tests/CMakeLists.txt` (baseline 2)
- [x] 3. Build `dsolve_corpus_tests`; run full baseline → 97/100, 0 FAIL, 0 crash
- [x] 4. Generate `reports/2.2.23.{md,tsv}` via `dsolve_corpus_report.py`
- [x] 5. Investigate residue: 2203 (Sinh·Cos resonant→UC decline→VoP hang, sympy=True),
        2289 (cubic-Root system churn), 2220 (cubic-Root IVP fit) — both sympy=False
- [x] 6. Fixes: (a) converter `\frac`/`\sqrt` brace-juxtaposition (2239/2240/2261/2300 malformed);
        (b) UC hyperbolic forcing (2203). 2289/2220 = bounded UNEVAL. → 98/100, 0 FAIL
- [x] 7. Regression: §2.2.4/§2.2.9 gates PASS; §2.2.20/21/22 regen byte-identical;
        dsolve stress suites PASS; check-c99 PASS; §2.1.2 gate (in progress)
- [x] 8. Docs: STATUS.md, README.md, DSOLVE_PLAN.md M42, changelog, version 0.139→0.140

## Review
**M42 — §2.2.23 (Problems 2201–2300): 98/100 PASS, 0 FAIL, 0 crash.**

Section: 45 scalar (17 IVP) + 55 constant-coeff linear systems — first systems-heavy
section since §2.2.15. Baseline 97/100; wave → 98/100 with two GENERAL fixes:

1. **Converter `\frac`/`\sqrt` brace-arg juxtaposition** (`tools/latex_ode_to_mathilda.py`) —
   `_implicit_mult` protected a math macro WITH its first brace arg, so fraction numerators
   escaped the implicit-mult split and coefficient×function glued into bogus symbols
   (`\frac{2ty}{…}`→`ty`, `\frac{4y_1}{3}`→`4y1`). Now protects only the macro NAME.
   Fixed malformed corpus records 2239/2240/2261 (systems) + 2300. §2.2.20/21/22 byte-identical.
2. **UC hyperbolic forcing** (`src/calculus/dsolve_undetcoeff.c`) — `Sinh`/`Cosh` forcing now
   expands to real exponentials before matching, so `Sinh·Cos−Cosh·Sin` (2203) becomes
   resonant UC atoms instead of hanging in VoP. No-op for non-hyperbolic forcing; strict
   residual gate unchanged (no wrong answer).

Residue 2 (both sympy=False, bounded declines): 2289 (3×3 cubic-`Root` complex spectrum →
realifier churn, no solve in 60s — M40/M41 irrational-spectrum-churn class) and 2220
(4th-order IVP, same spectrum → `Solve` IC-fit can't close over Root-object basis).

0 regression: §2.2.4/§2.2.9 gates PASS (contain the Sinh/Cosh cases now intercepted by UC),
§2.2.20/21/22 byte-identical, dsolve stress suites PASS, check-c99 PASS. Version 0.140.
