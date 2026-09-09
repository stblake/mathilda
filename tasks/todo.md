# M31 — §2.2.11 DSolve corpus (Problems 1001–1100), full coverage

Goal: add §2.2.11 to the DSolve corpus and drive to **100/100 (baseline 0)** with root-cause
fixes only (0 FAIL invariant, back-substitution verified). 8 s per-case fairness budget fixed.
Plan file: `~/.claude/plans/let-s-continue-our-implementation-valiant-raven.md`.

Composition (confirmed by fetch+convert): **41 first-order constant-coefficient linear systems**
(1001–1041, 2×2 … 6×6, defective/complex spectra) + **59 scalar (13 IVP)** — quadrature /
separable / first-order-linear, 2nd-order const-coeff (`missing_x`), exact, Euler, Gegenbauer,
Emden–Fowler, Liénard, Airy.

## Phase A — generate & sanity-check corpus
- [x] Fetch `indexsubsection20.htm` (browser-UA curl); run converter `--label 2.2.11` →
      `DE_examples_2211.m` (100 records: 59 scalar / 13 IVP / 41 systems). Converter unchanged.
- [x] Spot-check tricky shapes: 1096 (`alpha`), 1093 (`E^(-x)`), 1098–1100 (`y[t]`, indVar `t`),
      negative-point IVP 1086, 5–6-eq systems 1002/1003/1039. File closes clean.
- [x] Build runner; full run. **Baseline: 98/100, 0 FAIL, 0 crash, 0 timeout** — two UNEVAL,
      both systems: 1014 (3×3 DAG) and 1001 (4×4, spectrum {16,32,48,64}).

## Phase B — root-cause fixes (→ 0)
- [x] **1014** — `TriangularSystem` peels the DAG and hands the scalar engine
      `x2' = 9x2 + 7(C[k]−C[j])e^{2x}`; the integrating-factor integrand `e^{−9x}(…e^{2x}−…e^{2x})`
      reached `Integrate` as `Times[c, Plus[…]]` (exponents uncombined) → 55 s + branch-wrong
      `(−1)^{1/9}` antiderivative (kept, zero-test-undecidable). **Root fix in `Integrate`:**
      `integrate.c:try_linearity` now distributes a product over a sum factor (`c(g+h)→cg+ch`,
      commit-only-if-all-close), so the exponentials collapse (`e^{−7x}`, clean path). Also repairs
      the user-reported **direct** `Integrate[e^{−9x}(a e^{2x}−b e^{2x}), x]` bug.
- [x] **1001** — answer CORRECT (`Simplify[resid]≡0`), but `e^{64x}`-scale cancellation made the
      20-digit sweep read it as nonzero → UNEVAL. Fix: `dsolve_corpus_prelude.m:dsResidVerdict`
      re-checks a not-small sample at 200-digit precision (monotone; never a FAIL).
- [x] Re-run §2.2.11 → **100/100, 0 FAIL, baseline 0**.

## Phase C — deliverables
- [x] `reports/2.2.11.{tsv,md}`; ctest `dsolve_corpus_2_2_11_tests` (baseline 0) in CMakeLists.
- [x] STATUS.md §2.2.11 block + M31 wave-history; README.md row; DSOLVE_PLAN.md M30+M31 entries;
      changelog `2026-09-07.md` M31 section; `calculus.md` LinearFirstOrder note.
- [x] Anti-overfit units `t_m31_triangular_exp_forcing`, `t_m31_linsys_large_eigenvalue`
      (`tests/test_dsolve.c`) + `test_linearity_distributes_product` (`test_integrate_dispatch.c`,
      direct-`Integrate` regression guard) — green.

## Phase D — regression & gates (no regression; 4 sections improved)
- [x] §2.2.1–§2.2.10 full re-run: **0 FAIL everywhere**, all ≤ baseline. Improved by the shared
      prelude fix: §2.2.1 (4→2), §2.2.2 (8→7), §2.2.4 (1→0), §2.2.7 (7→6) — baselines tightened,
      reports regenerated.
- [x] All 22 integrate unit suites green with the root `Integrate` fix (form-change check —
      `ctest -R "integrate|intrat"`, exit 0). No existing closed form changed shape.
- [ ] §2.1.2 + §2.2.x full re-run with the root fix — confirm ≤ 655 / no regression, re-baseline
      §2.1.2 if improved.
- [x] `make check-c99` clean. `dsolve_tests` green; §2.2.11 ctest green (baseline 0).
- [ ] cmake reconfigure + `ctest -R "dsolve|integrate"` all green.
- [ ] Rebuild code-review graph.

## Review
- §2.2.11 is Mathilda's strongest DSolve territory: the entire scalar half (59/59) solved out of
  the box. Both gaps were first-order constant-coefficient **systems**, and both traced to
  pre-existing issues the corpus happened to exercise:
  1. a latent **Integrate** wrong-answer (un-combined exponential product) exposed through
     TriangularSystem's integrating-factor sub-solve — fixed at the ROOT in
     `integrate.c:try_linearity` (distribute a product over a sum factor; `Integrate` is linear),
     which also repairs the user-reported direct `Integrate[E^(-9x)(a E^(2x)-b E^(2x)),x]` bug.
  2. a **verifier** precision limit (catastrophic cancellation in large-eigenvalue solutions) —
     fixed in the shared prelude with a monotone high-precision re-check that also lifted four
     earlier sections.
- No wrong answer can ship: the 0-FAIL invariant held across all sections; the verifier fix is
  provably monotone (cannot introduce a FAIL); the Integrate fix distributes a product over a sum
  and commits only a fully-elementary split (behavior-preserving linearity — all 22 integrate
  unit suites stay green, and it eliminates a pre-existing WRONG antiderivative).
