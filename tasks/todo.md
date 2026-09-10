# M35 — DSolve §2.2.15 corpus (Problems 1401–1500) + piecewise/step forcing

## 1. Converter (`tools/latex_ode_to_mathilda.py`)
- [x] `replace_cases` pre-pass: `\left\{…cases…\right.` → `Piecewise[{{v,(c)},…},0]` (protect before split; detect symbols on UNprotected tex so `t` isn't hidden)
- [x] `convert_side`: `\le`/`\leq`→`<=`, `\ge`/`\geq`→`>=`, `\infty`→`Infinity`; drop `±Infinity` bounds in clause conditions
- [x] `Heaviside` → `UnitStep`
- [x] Regenerate `DSolve_test_status/DE_examples_2215.m`; all 100 records parse

## 2. Baseline measurement
- [x] Registered `dsolve_corpus_2_2_15_tests`; baseline 98/100, 0 FAIL, 2 UNEVAL (1463,1469)
- [x] `reports/2.2.15.{md,tsv}` generated; STATUS.md updated
- [x] Discovered 1492-1500 FALSE-passed at baseline (inert Integrate[UnitStep…], UNK→trusted)

## 3. New method `DSolve`PiecewiseForcing` (`src/calculus/dsolve_piecewise.c`)
- [x] Detect linear scalar IVP with piecewise/step forcing; breakpoints via UnitStep root + Piecewise cond operands
- [x] Interval continuation (recurse DSolve per piece; continuity handoff; pw_resolve per interval)
- [x] Assemble `Piecewise[...]`; `pw_num_ok` internal guard (residual per interval + each IC)
- [x] Bounded/re-entry (TimeConstrained + 8s deadline + decline memo)
- [x] Cascade slot before `dsolve_undetcoeff_try`; registered in `dsolve_init`
- [x] Shared fix: `ds_residual_is_distributional` +UnitStep/Piecewise (zero_test mis-rejects piecewise residuals; e.g. 1497)

## 4. Hard residues
- [x] 1463 (4th-order transcendental), 1469 (3rd-order variable-coeff): documented bounded declines (both also ✗ in SymPy)

## 5. Wiring & dashboard
- [x] ctest baseline = 2 (the true residue)
- [x] STATUS.md §2.2.15 block; README contents rows (2214 + 2215)
- [x] `t_m35_piecewise_forcing` in `test_dsolve.c` (+ CMake COMMON_SRC entry)
- [x] DSOLVE_PLAN.md M35 entry; version bump 0.132→0.133
- [x] Changelog note (docs/spec/changelog/2026-09-07.md)

## 6. Verification
- [x] `make -j` clean (no warnings) + `make check-c99` green
- [x] **All 15 §2.2.1–§2.2.15 corpus gates PASS (100%, 0 failed)** — no regression
- [x] valgrind: solve works, 13.4KB leak == baseline plain-IVP (inherited engine leak; dsolve_piecewise adds none)
- [x] Confirmed pre-existing: dsolve_tests SIGALRM at t_rischnorman (Abel hang) reproduces on clean M34 tree
- [x] rebuilt code-review graph

## Review

**M35 complete.** §2.2.15 (Boyce & DiPrima 1401–1500) added to the corpus at 98/100, 0 FAIL.
The headline: the 9 step/piecewise-forced IVPs (1492–1500) that FALSE-passed at baseline
(inert `Integrate[UnitStep…]`, scored UNK→trusted) now GENUINELY solve, via the new
`DSolve\`PiecewiseForcing` (interval continuation → verified `Piecewise`, Mathematica's form),
directly verified (residual ~0 in every interval + each IC). Two shared substrate fixes:
`ds_residual_is_distributional` accepts UnitStep/Piecewise residuals (zero_test mis-rejects them),
and the converter learned the `cases`/`Heaviside` LaTeX. Residue 2 (1463, 1469) are research-grade
higher-order variable-coefficient equations SymPy also fails. No regression (corpus gates hold;
valgrind flat vs baseline; the dsolve_tests SIGALRM is a pre-existing Abel-hang, confirmed on the
clean M34 tree). Version 0.133.
