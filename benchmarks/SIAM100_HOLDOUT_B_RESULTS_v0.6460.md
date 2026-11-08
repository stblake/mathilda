# SIAM 100-digit: sealed holdout B, the single run at v0.6460

Holdout B is experiments 226–255: 30 families, 10 problems each, at two tiers (machine precision and 50 digits), so 600 cases. It was written by an agent that had seen neither holdout A nor any fix, and the problems differ from holdout A in kind where they could. Every truth value has two independent sources (mpmath and PARI, or a closed form), recorded in `benchmarks/siam100_refs/b*`.

The other four systems' baselines were recorded first, on a quiet machine. Mathilda then ran **once**, at v0.6460 (commit ab0e54c6e), on 2026-11-09. No fix was made against holdout B before or since; from now on it is training data. Measuring generalisation again needs a fresh holdout C.

- Full tables: `siam100_results/holdoutB_v0.6460.md` (all systems) and `siam100_results/holdoutB_baseline.md` (the other systems before Mathilda ran).
- Raw JSONs: `siam100_results/raw_v0.6460.tar.xz`.

The scoring is as in `SIAM100.md`. A case is "on target" when Mathilda has at least as many correct digits as the best other system (counted with the book's truncation rule) **and** takes less than 1.5× the time of the fastest correct system.

## Headline

| set | Mathilda version | in or out of sample | on target |
|---|---|---|---|
| holdout A, before any fix | v0.6377 | out | 262 / 600 (43.7 %) |
| **holdout B** | **v0.6460** | **out** | **395 / 600 (65.8 %)** |
| holdout A, after the fixes | v0.6460 | in (training) | 506 / 600 (84.3 %) |
| original 30 problems (76 cases) | v0.6460 | in | 69 / 76 |

The campaign's fixes raised the out-of-sample on-target rate by 22 percentage points, from 43.7 % to 65.8 %. That is 133 of the 244 cases gained on holdout A. The other 111 cases, the gap between 84.3 % and 65.8 %, are where a fix covered the family members holdout A happened to contain and not the variants in holdout B.

Classes per system on holdout B:

| system | correct | inexact | wrong | none |
|---|---|---|---|---|
| **Mathilda 0.6460** | **395** | 79 | 45 | 81 |
| Python (scipy/mpmath) | 318 | 69 | 44 | 169 |
| WL | 275 | 213 | 61 | 51 |
| Maxima | 185 | 26 | 11 | 378 |
| FriCAS | 108 | 16 | 3 | 473 |

Mathilda has the most correct answers of the five systems: 395 against 318 for the next best, Python. Its wrong count, 45, is level with Python's 44 and below WL's 61. 20 of the 45 come from a single defect in D17 (below). Without it, Mathilda's wrong count would be 25.

## Per family: on target out of 20

| family | A before fixes (v0.6377) | A after fixes (v0.6460) | **B (v0.6460)** |
|---|---|---|---|
| C1 oscillatory tails | 10 | 20 | 17 |
| C2 billiards | 12 | 13 | 13 |
| C3 matrix norms | 10 | 18 | 13 |
| C4 global minima | 13 | 20 | **20** |
| C5 discrete minimax | 7 | 19 | 14 |
| C6 lattice walks | 8 | 13 | 6 |
| C7 sparse solves | 7 | 20 | 18 |
| C8 heat times | 2 | 20 | 19 |
| C9 parametric optima | 9 | 19 | 14 |
| C10 exit probabilities | 12 | 20 | **20** |
| D1 digit-restricted series | 17 | 17 | 16 |
| D2 digit-length series | 16 | 16 | 6 |
| D3 minimal energy | 8 | 19 | 19 |
| D4 series roots | 12 | 19 | 18 |
| D5 matrix-polynomial norms | 17 | 17 | 14 |
| D6 oscillatory towers | 8 | 14 | 9 |
| D7 region areas | 3 | 19 | 18 |
| D8 rational ℝⁿ integrals | 17 | 17 | 7 |
| D9 lattice walks | 10 | 16 | 18 |
| D10 random recurrences | 1 | 6 | 0 |
| D11 vector shooting | 1 | 19 | 17 |
| D12 event-time roots | 1 | 18 | 15 |
| D13 initial-value thresholds | 14 | 15 | 13 |
| D14 FEM parameter | 6 | 16 | 10 |
| D15 blow-up times | 8 | 11 | 13 |
| D16 heat functionals | 12 | 20 | 17 |
| D17 nodal lines | 14 | 20 | **0** |
| D18 Dirichlet eigenvalues | 2 | 20 | 14 |
| D19 design | 4 | 15 | 11 |
| D20 functional roots | 1 | 10 | 6 |

**The fixes generalised** where they were made in the generic head (NIntegrate, FindRoot, NDSolve events, NMinimize/FindMinimum, sparse LinearSolve, NSum):
- C4, C7, C8, C10, D3, D4, D7, D11 and D12 are within 0–3 cases of their holdout-A result;
- D9 and D15 are better than on holdout A.

**The fixes did not generalise in these families:**

- **D17 (20 → 0): every case is wrong, with the same value.** Mathilda returns `{0, 0}` for the nodal-line length on every domain, in about 0.3–0.4 s. One systematic defect, not 20 separate ones. Not yet diagnosed: either the eigenfunction step or the nodal-set measurement finds no nodal segment on these domains.
- **D8 (17 → 7):** the rational integrands over ℝⁿ time out, 10 of them, on the multi-dimensional infinite-range `NIntegrate` path.
- **D2 (16 → 6):** `NSum` stops with `NSum::accgl` on these digit-length weights. Mathilda answers in only one of the ten 50-digit cases, and two answers are wrong.
- **D14 (16 → 10):** 8 cases fail with `NDSolve::femfail` ("not linear in the unknown") or `NDSolve::femspec`. FEM accepts only the linear operator forms holdout A happened to use.
- **D10 (6 → 0):** every answer is correct, but 3–13× slower than the fastest system, and a few are a digit or two short. It is the same `FoldList` speed gap as on holdout A.
- **C6 (13 → 6):** the 50-digit 2-D and 3-D integrals time out. This was already open on holdout A.
- **D18 and D19 (6 and 8 unevaluated):**
  - D18: `NDEigenvalues` on region forms holdout A did not contain;
  - D19: the constrained design problems.
- **D20 (4 unevaluated), C3 and C5 (correct but slower):**
  - C3: SVD at 1.9–4.2× the fastest system;
  - C5: minimax at 10–20×.

## What this means

The campaign improved Mathilda well beyond the problems it was trained on. That is the honest claim, and the sealed holdout B is the evidence. Mathilda now has the most correct answers of the five systems on unseen problems. It is not yet "at least as good as every other system on every case": it is so in two thirds of unseen cases, and in 84 % of the cases it has seen. The list above is the next work queue. Fixing any of it turns holdout B into training data, so the next measurement needs a holdout C.
