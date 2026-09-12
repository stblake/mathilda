# M39 — §2.2.19 corpus (Problems 1801–1900) + push to 100/100

Baseline measured this session: **94/100, 0 FAIL, 6 UNEVAL**. Target: 100/100,
0 FAIL, 0 regression on prior sections.

## Step 1 — Land corpus + harness (no solver change)
- [ ] Generate `DSolve_test_status/DE_examples_2219.m` (--label 2.2.19, --url indexsubsection28)
- [ ] Verify §2.2.1–§2.2.18 regenerate byte-for-byte (no converter change)
- [ ] Register `dsolve_corpus_2_2_19_tests` in `tests/CMakeLists.txt` (baseline 6)
- [ ] Generate `reports/2.2.19.{tsv,md}`; add §2.2.19 block to STATUS.md + README row

## Step 2 — Cluster A: change-of-variable + VoP (1817, 1822, 1823)
- [ ] Power substitution `t=x^k` in `dsolve_changevar.c` (→ Bessel; 1817 homog)
- [ ] Carry forcing through changevar → VoP + numeric verify (1822, 1817)
- [ ] Bound/robustify `dsolve_variation_of_parameters` (1823 timeout)

## Step 3 — Cluster B: IVP fit + form normalization (1836, 1876)
- [ ] Normalize `E^(Σ cᵢ Log[gᵢ])` → `∏ gᵢ^cᵢ` (1876 branch artifacts)
- [ ] IVP fitter robustness for special-function bodies (1836 Ei-form)

## Step 4 — Cluster C: exact 2nd-order (1900)
- [ ] `dsolve_exactode.c`: first-integral / implicit form instead of series fall

## Step 5 — Re-measure, lower baseline, document
- [ ] Re-run §2.2.19 + all prior sections (0 regression); regen reports
- [ ] Update STATUS.md counts + wave history; lower baselines in CMakeLists
- [ ] Anti-overfit unit tests in `tests/test_dsolve.c`
- [ ] `make check-c99`; DSolve ctest + stress green; valgrind spot-check
- [ ] Docs: M39 in DSOLVE_PLAN.md; changelog 2026-09-07; version 0.136→0.137
- [ ] Rebuild code-review graph

## Review

**Outcome: §2.2.19 94/100 → 95/100, 0 FAIL, 0 regression (M39).** Version 0.136 → 0.137.

- **Step 1 (corpus)**: `DE_examples_2219.m` (100 scalar, 27 IVP, 0 systems) landed via the converter
  (no converter change); `dsolve_corpus_2_2_19_tests` registered, baseline lowered 6 → 5;
  `reports/2.2.19.{tsv,md}`, STATUS.md §2.2.19 block, README row all added.
- **Step 2 (Cluster A)**: new method **`DSolve\`VariationOfParameters`** (`dsolve_nonhomog_vop.c`) —
  the transcendental-coefficient nonhomogeneous backstop (solve homogeneous → normalise fundamental
  set with `PowerExpand[Simplify[·]]` → VoP → numeric-verify). Closes **1822**. Anti-overfit
  `t_m39_nonhomog_vop`.
- **Steps 3–4 not landed (documented residue)**: 1876 (Kovacic `E^(nLog)` IC-fit branch artifacts —
  a form-normalisation gap; the fix touches `assemble_general`, regression-sensitive), 1836 (Solve
  can't fit an `ExpIntegralEi`-constant IVP system), 1823 (Kovacic constant-`r` churn before the
  gated backstop is reached), 1817 (non-elementary Struve particular), 1900 (exact → non-elementary /
  complex-singular ₂F₁). All bounded declines, no wrong answers.

**Key lessons (also in harness memory):**
- A DSolve method that recursively calls `DSolve`/`Simplify` must NOT wrap the sub-call in
  `TimeConstrained` — nested TimeConstrained aborts the whole subtree (`$Aborted`). Bound with a
  `time()`-deadline instead.
- A "solve-homogeneous-then-VoP" backstop belongs AFTER the closed-form 2nd-order methods (pure
  backstop → 0-regression by construction) and must be tightly gated (transcendental coefficients)
  or it slows the whole corpus with redundant recursive re-solves; placing it BEFORE Kovacic
  regressed 1828/1835.
- Radical-of-trig/power basis cleanup order is **Simplify then PowerExpand** (Simplify canonicalises
  `1+Tan²→Sec²`; PowerExpand reduces `1/Sqrt[Sec²x]→Cos x`); Simplify AFTER PowerExpand reverts it.
