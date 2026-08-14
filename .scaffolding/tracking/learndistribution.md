---
schema_version: 2
title: Machine learning family 4 — LearnDistribution
slug: learndistribution
status: implementation
source: direct-user-request
owner: Michael Sollami
issue: pending
pull_request: pending
started: 2026-08-14
last_updated: 2026-08-14
blocked_by: none
goal_lock:
  status: active
  stamped: 2026-08-14 02:35
  scope:
    - "src/ml/**"
    - "src/random.c"
    - "src/random.h"
    - "src/eval.c"
    - "src/core.c"
    - "src/sym_names.c"
    - "src/sym_names.h"
    - "src/print.c"
    - "makefile"
    - "tests/**"
    - "docs/spec/**"
  success_criteria:
    - "AC-1 — RandomVariate is reproducible under SeedRandom, including across an odd number of prior draws."
    - "AC-2 — Sample moments match the distribution (mean and standard deviation within tolerance at n = 20000)."
    - "AC-3 — PDF agrees with the closed form computed symbolically."
    - "AC-4 — Invalid parameters decline rather than producing NaN."
    - "AC-5 — LearnDistribution fits Multinormal and GaussianMixture and exposes a density."
---

# Machine learning family 4 — LearnDistribution

## Feature Definition

- One-line goal: `Implement LearnDistribution and the distribution/sampling substrate it needs.`
- Problem: `None of LearnDistribution, RandomVariate, PDF, NormalDistribution, MultinormalDistribution, SmoothKernelDistribution or ContingencyTable existed. More pointedly there was NO GAUSSIAN DEVIATE anywhere in the tree -- RandomReal is uniform and the xoshiro stream under it produces uniforms -- so no distribution-based method could sample at all.`
- Requested by: `direct user request, continuing the ML families`

## Current State Study

- `random_uniform_01` (src/random.c:609) was `static`; it is the correct base for
  deviates because it draws from the USER-VISIBLE stream, which is what makes
  SeedRandom apply. `random_internal_int_range` is explicitly the wrong path — it is
  reserved for internal decision procedures that must NOT perturb that stream.
- `ml_gmm_fit` / `ml_gmm_logpdf` already exist in `src/ml/gmm.h`, extracted in
  iteration 9 for exactly this family.
- The trained-model representation (src/ml/predict.h) is available for a fitted
  distribution, but a SPECIFIED distribution is a different thing and must not reuse
  it — see Decisions.

## Scope

### In Scope

- `RandomVariate`, `PDF`, `NormalDistribution`, `UniformDistribution` (done, iteration 16)
- `LearnDistribution` with `"Multinormal"` and `"GaussianMixture"`
- `MultinormalDistribution`, and `SmoothKernelDistribution` if it fits

### Out Of Scope

- The full Wolfram distribution zoo (Poisson, Binomial, Beta, ...) — this is the ML
  family, not a statistics library
- `ContingencyTable` unless it lands cheaply after the above

## Success Criteria

- `AC-1` — reproducible under SeedRandom, including after an odd number of draws
- `AC-2` — sample moments match
- `AC-3` — PDF matches the closed form
- `AC-4` — invalid parameters decline
- `AC-5` — LearnDistribution fits and exposes a density

## Tests

- `tests/test_ml_dist.c`, registered with BOTH `add_executable` and `add_test`.
  `random_tests` is in this family\'s regression set because `random.c` changed.

## Task List

- [x] `Expose random_uniform_01; Box-Muller Gaussian deviate; RandomVariate; PDF; Normal and Uniform distribution objects` | `independent` | `done (iteration 16)`
- [x] `LearnDistribution with Method -> "Multinormal"` | `depends on: distributions` | `done (iteration 17)`
- [x] `LearnDistribution with Method -> "GaussianMixture"` | `depends on: distributions` | `done (iteration 18) -- the payoff for the iteration-9 extraction`
- [x] `SmoothKernelDistribution / KernelDensityEstimation` | `depends on: distributions` | `done (iteration 19)`
- [x] `ContingencyTable` | `independent` | `OUT OF SCOPE -- recorded, not implemented. It is a contingency-table builder over CATEGORICAL data, which needs the categorical encoding that family 5 (Classify) must design anyway. Building it first would mean inventing that encoding twice. Deferred to family 5 as a follow-up rather than done badly here.`

## Checkpoints

- [x] start | completed: `2026-08-14 02:35`
- [x] spec / plan created | completed: `2026-08-14 02:35`
- [x] implementation started | completed: `2026-08-14 02:35`
- [ ] implementation complete | completed: `pending`
- [ ] critic pass | completed: `pending`
- [ ] risk-register reviewed | completed: `pending`
- [ ] feature validated | completed: `pending`
- [ ] PR created | completed: `pending`
- [ ] closeout complete | completed: `pending`

## Decisions

- `A SPECIFIED distribution is NOT the trained-model representation, even though the two look alike. A distribution is written by the user, so its parameters ARE the information and it prints in full; a fitted model derives its parameters and elides them. Same mechanism, opposite convention on visibility, deliberately not unified -- a test row pins the full printing so a later change cannot merge them by accident.`
- `Sampling goes through random_uniform_01 (the user-visible stream) rather than a private generator, so SeedRandom applies. A private generator would leave reproducibility HALF working -- honoured by RandomReal, ignored by RandomVariate -- which is worse than not working.`
- `Box-Muller in the POLAR form: no sin/cos, and rejecting the corners of the square is cheaper than two transcendental calls.`
- `SeedRandom clears the cached spare deviate. Without this the first draw after reseeding comes from the previous stream, and the bug only appears after an ODD number of prior draws, so an even-count test would hide it. Both parities are asserted.`

## Risks And Unknowns

- `A sampler is hard to test because its output is meant to be unpredictable. Mitigated by asserting the three deterministic things: reproducibility, moments at large n, and the closed-form PDF.`
- `Moment tolerances could be flaky. 20000 samples puts the standard error of the mean at sigma/141, so the 0.1 tolerance on sigma = 2 is about seven standard errors -- loose enough not to flake, tight enough to catch a wrong variance.`

## Risk Register

| ID | Category | Likelihood | Impact | Mitigation | Residual | Owner |
|---|---|---|---|---|---|---|
| RISK-1 | Half-working reproducibility | high | high | A cached Gaussian deviate surviving SeedRandom makes RandomVariate irreproducible while RandomReal stays reproducible — the failure is invisible in the common test. Cache cleared in builtin_seedrandom; both draw parities asserted. | Closed. | Michael Sollami |
| RISK-2 | Silent NaN propagation | med | med | A non-positive sigma would produce NaNs through an entire sample, surfacing much later as a strange plot. Invalid parameters decline instead. | Closed. | Michael Sollami |
| RISK-3 | Flaky statistical tests | med | low | Tolerances set at ~7 standard errors and the seed fixed, so the test is deterministic in practice. | Accepted. | Michael Sollami |

## Test Results

- Command: `cd tests/build && ./ml_dist_tests` plus ml_pca, ml_predict, random, eval, print, findclusters, list
- Outcome: `all pass (7 ml_dist tests); check-c99 / check-packed-aware / check-array-exactness OK`

## Proofs Of Completion

- Completion timestamp: `pending`

## Tech Debt Review

- Potential tech debt introduced:
  - `The Box-Muller cache is per-process static state. Correct here, but it would need care if the evaluator ever became multi-threaded.`
  - `RESOLVED (iteration 17): the mu/chol/dim fields are now used by the Multinormal path.`
- Existing tech debt noticed:
  - `PDF supports only the two distributions implemented here; Wolfram PDF is generic over its whole distribution zoo.`
- Mitigations taken:
  - `Both draw parities asserted for the cache reset; invalid parameters decline.`
- Follow-up needed:
  - `LearnDistribution Multinormal and GaussianMixture; SmoothKernelDistribution; ContingencyTable.`

## Activity Log

- `2026-08-14 03:20` Iteration 19: SmoothKernelDistribution, and FAMILY 4 IS CLOSED. Bandwidth is the multivariate normal-reference rule, which reduces exactly to Silverman in 1-D; product (diagonal) kernel because a full-covariance kernel needs a bandwidth MATRIX, a much harder estimation problem than the default rule solves. The verification is an EXACT identity rather than an approximation: a KDE is the empirical distribution convolved with the kernel, so Var = ML sample variance + h^2 and Mean = sample mean, which pins the bandwidth rule, the kernel variance and the normalisation in one assertion. It first failed at 1.6e-5, and refining the integration grid did NOT shrink the error -- that constancy is what identified it as a wrong CONSTANT in my test (the rounded 1.0592 instead of the exact (4/3)^(1/5) = 1.059224) rather than quadrature error. With the exact constant it holds to 3.8e-15. ContingencyTable judged OUT OF SCOPE here: it is categorical, and the encoding it needs is exactly what family 5 must design, so building it now would mean inventing that twice.
- `2026-08-14 03:05` Iteration 18: LearnDistribution with GaussianMixture -- the payoff for extracting ml_gmm_fit in iteration 9, and it was indeed small. Variance floor for a standalone fit: the SQUARED MEDIAN NEAREST-NEIGHBOUR DISTANCE, the standalone analogue of fc_gmm_ndim's median spanning-tree edge weight; median not mean so one tight pair cannot drag it toward zero, which matters because a mixture's likelihood is unbounded above and a merely-small floor lets BIC buy near-singular spikes. The k=1 cross-check against Multinormal turned out NOT to be exact equality, and pinning down WHY made the test stronger than the one I had planned: EM uses the ML (n) divisor, Multinormal the unbiased (n-1) one, and the mixture adds the ridge, so cov_mix = ((n-1)/n) cov_mn + floor -- measured 1.41 = (40/41)(1.435) + 0.01 exactly, with the floor 0.01 being the squared 0.1 spacing. A fuzzy 'ratio about 1.009' assertion would have hidden a wrong divisor OR a wrong floor; this pins both. Chol and logdet are recomputed on read rather than stored, so the same fact does not live in two places.
- `2026-08-14 02:50` Iteration 17: LearnDistribution with Multinormal. Promoted ml_chol and ml_mahalanobis out of gmm.c into src/ml/mlutil -- second real consumer, the same rule that promoted the list builders and that declined to extract find_clusters.c's metric layer. Verified the GMM still passes after the move before going further. The independent cross-check is strong here: a 1-D fit reaches its density via Cholesky + Mahalanobis while PDF[NormalDistribution] uses the scalar closed form, sharing no code, and they agree exactly including 3.5 sigma into the tail where a wrong log-determinant would show only as a small relative error. Added a normalisation test too -- the density integrates to 1.0 over +/-6 sigma -- because an agreement test alone would pass two densities sharing a normalisation error. One asymmetry worth noting in the code: for a multinormal a flat list is ONE point (its argument is a list), the opposite of the scalar case where a list is many points.
- `2026-08-14 02:35` Created, and iteration 16 landed: random_uniform_01 exposed, a Box-Muller Gaussian deviate (the tree had none), RandomVariate, PDF, and the Normal/Uniform distribution objects. The non-obvious bug found and fixed by design rather than by accident: SeedRandom must clear the cached spare deviate, or reproducibility half-works -- honoured by RandomReal, silently broken for RandomVariate -- and it only shows after an ODD number of prior draws, so both parities are asserted. Cross-checked PDF against the closed form computed by the existing symbolic machinery, which shares no code with the C density.

## Reflection

- What went well: `pending`
- What went wrong: `pending`
- Gaps to close: `pending`
- Skill / AGENTS.md updates to propose: `pending`

## Follow-Up

- `Family 5: Classify + NaiveBayes, NearestNeighbors, LogisticRegression, DecisionTree, RandomForest. Needs feature types and categorical encoding, which is the real work.`

## Team Addendum

- `none`
