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
- [ ] `LearnDistribution with Method -> "Multinormal"` | `depends on: distributions` | `pending`
- [ ] `LearnDistribution with Method -> "GaussianMixture" (wire up ml_gmm_fit/ml_gmm_logpdf)` | `depends on: distributions` | `pending`
- [ ] `SmoothKernelDistribution / KernelDensityEstimation` | `depends on: distributions` | `pending`
- [ ] `ContingencyTable` | `independent` | `pending`

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
  - `MlDist carries unused mu/chol/dim fields for a Multinormal that is not implemented yet -- they land with the next piece.`
- Existing tech debt noticed:
  - `PDF supports only the two distributions implemented here; Wolfram PDF is generic over its whole distribution zoo.`
- Mitigations taken:
  - `Both draw parities asserted for the cache reset; invalid parameters decline.`
- Follow-up needed:
  - `LearnDistribution Multinormal and GaussianMixture; SmoothKernelDistribution; ContingencyTable.`

## Activity Log

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
