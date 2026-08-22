---
schema_version: 2
title: Machine learning family 2 — dimensionality reduction
slug: dimensionreduce
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
  stamped: 2026-08-14 01:30
  scope:
    - "src/ml/**"
    - "src/list/**"
    - "src/linalg/**"
    - "src/core.c"
    - "src/sym_names.c"
    - "makefile"
    - "tests/**"
    - "docs/spec/**"
    - "tools/*.py"
  success_criteria:
    - "AC-1 — Standardize agrees with (x - Mean[x])/StandardDeviation[x] written by hand, including the n-1 divisor."
    - "AC-2 — PrincipalComponents preserves total variance (the rotation is orthogonal) and orders components by decreasing variance."
    - "AC-3 — Rank-deficient input puts zero variance in the trailing components."
    - "AC-4 — Method -> Correlation differs from Covariance where the columns have different scales."
    - "AC-5 — Every new builtin returns head List, not a visible NDArray, matching Inverse/Dot/LinearSolve."
    - "AC-6 — Eigenvector signs are canonical, so results do not depend on whether LAPACK was linked."
---

# Machine learning family 2 — dimensionality reduction

## Feature Definition

- One-line goal: `Implement DimensionReduce and its methods, plus the PrincipalComponents and Standardize primitives they rest on.`
- Problem: `Wolfram's Machine Learning Methods guide lists DimensionReduction with ten methods. Mathilda had NONE of DimensionReduce, PrincipalComponents, Standardize or FeatureExtract -- verified by grepping src/ for the quoted names. SVD and symmetric eigen do exist at machine precision, LAPACK-backed, so the linear-algebra substrate is there and the gap is the statistical layer on top.`
- Requested by: `direct user request, continuing the ML families after FindClusters`

## Current State Study

- `numarray.h` EXPORTS the machine bridge (`na_load_matrix`, `na_build_matrix`) — the
  right loader, already public, so no new one was written.
- `na_build_matrix` returns a **visible NDArray** (head `NDArray`). `Inverse`, `Dot`
  and `LinearSolve` return head `List`. Using the bridge's builder for a new
  user-facing result therefore breaks `result === {{...}}`, which is why these
  builtins construct plain Lists and let the evaluator's packing gate decide.
- `mat_lapack_dsyev` is declared unconditionally; without LAPACK a stub returns
  nonzero. So no `#ifdef` is needed at a call site — check the status and fall back.
- Attributes are set as `symtab_get_def(name)->attributes |= ATTR_PROTECTED`; there is
  no `symtab_set_attributes`.
- `Expr`'s string payload is `e->data.string` (a `char*`), not `e->data.string.value`.

## Scope

### In Scope

- `Standardize`, `PrincipalComponents` (done, iteration 10)
- `DimensionReduce` with a method dispatch, plus `MultidimensionalScaling` and
  `LatentSemanticAnalysis`
- The reusable buffer-level kernels in `src/ml/pca.{c,h}`

### Out Of Scope

- `FeatureExtract` / feature-type inference — belongs with Classify (family 5)
- The remaining seven Wolfram DimensionReduction methods (autoencoders, t-SNE, UMAP)
- A trained-reducer object that can be applied to new data — needs the model
  representation being deferred to family 3

## Success Criteria

- `AC-1` — `Standardize` matches its definition including the n-1 divisor
- `AC-2` — PCA preserves total variance and orders components by decreasing variance
- `AC-3` — rank-deficient input gives zero trailing components
- `AC-4` — `Correlation` differs from `Covariance` on differently-scaled columns
- `AC-5` — results are head `List`, not visible `NDArray`
- `AC-6` — eigenvector signs are canonical across backends

## Tests

- `tests/test_ml_pca.c`, registered with BOTH `add_executable` and `add_test`.
  Assertions are properties any correct implementation must satisfy (definitional
  agreement, variance preservation, rank deficiency, option-changes-answer) rather
  than numbers this implementation happens to produce — there is no prior behaviour to
  pin, since these are new builtins rather than ports.

## Task List

- [x] `src/ml/pca.{c,h} kernels: column mean/sd, standardise, symmetric eigen desc, PCA` | `independent` | `done (iteration 10)`
- [x] `Standardize + PrincipalComponents builtins, registered with attributes and docstrings` | `depends on: kernels` | `done (iteration 10)`
- [x] `DimensionReduce with a method dispatch (PCA first)` | `depends on: kernels` | `done (iteration 11)`
- [x] `MultidimensionalScaling` | `depends on: symmetric eigen` | `done (iteration 11, as a DimensionReduce method)`
- [x] `LatentSemanticAnalysis` | `depends on: SVD` | `done (iteration 11, as a DimensionReduce method -- the Gram-matrix eigendecomposition IS a truncated SVD, so no separate SVD path was needed)`
- [x] `A DimensionReducerFunction applicable to new data` | `depends on: the family-3 model representation` | `done (iteration 14) -- DimensionReduction[data, k]; family 2 is now CLOSED`

## Checkpoints

- [x] start | completed: `2026-08-14 01:30`
- [x] spec / plan created | completed: `2026-08-14 01:30`
- [x] implementation started | completed: `2026-08-14 01:30`
- [ ] implementation complete | completed: `pending`
- [ ] critic pass | completed: `pending`
- [ ] risk-register reviewed | completed: `pending`
- [ ] feature validated | completed: `pending`
- [ ] PR created | completed: `pending`
- [ ] closeout complete | completed: `pending`

## Decisions

- `Kernels are buffer-level (row-major n x dim double) in src/ml/, following the precedent set by src/ml/gmm.c. Column statistics are wanted by Standardize, PCA, a future Classify's feature scaling and LearnDistribution's Multinormal; the symmetric eigendecomposition is wanted by PCA, MDS and LSA. Writing them against Expr would put all of that behind a builtin.`
- `Results are plain Lists, not na_build_matrix's visible NDArray, so that === against a literal behaves the way it does for Inverse and Dot. This was found by comparing surfaces, not assumed.`
- `Eigenvector signs are canonicalised (largest-magnitude component positive). Without this the same input gives sign-flipped output depending on whether LAPACK was linked, and nothing could be pinned in a test.`
- `Jacobi is the no-LAPACK fallback and is written in-house, for the reason the Cholesky in gmm.c is: the in-tree symmetric solvers are static to their translation units, dim is a feature count, and this runs once per call.`

## Risks And Unknowns

- `PCA on the correlation matrix silently changes what "variance explained" means; a caller comparing eigenvalues across the two methods is comparing different quantities.`
- `MDS stress and PCA variance ratios are DERIVED quantities, and the Spectral port showed that a statistic correct on raw data can be wrong on a derived one (median embedding jump -> 0). Expect the same class of trap.`

## Risk Register

| ID | Category | Likelihood | Impact | Mitigation | Residual | Owner |
|---|---|---|---|---|---|---|
| RISK-1 | Surface inconsistency | high | med | A new builtin using the machine bridge's builder returns a visible NDArray and compares False against a literal list. Caught by comparing against Inverse/Dot; a test row now pins head == List. | Closed for these two builtins; the trap remains for any future one. | Michael Sollami |
| RISK-2 | Build-dependent output | med | high | Eigenvector signs differ between LAPACK and Jacobi. Canonicalised in ml_sym_eigen_desc so results are reproducible and testable. | Closed. | Michael Sollami |

## Test Results

- Command: `cd tests/build && ./ml_pca_tests` plus the four FindClusters suites and linalg_tests
- Outcome: `all pass (9 ml_pca tests); make check-c99 / check-packed-aware / check-array-exactness OK`

## Proofs Of Completion

- Completion timestamp: `pending`

## Tech Debt Review

- Potential tech debt introduced:
  - `ml_list_of_reals / ml_list_matrix are local to pca.c. Every future src/ml builtin needs the same thing, so the second consumer should promote them rather than copy them.`
- Existing tech debt noticed:
  - `na_build_matrix returning a visible NDArray is a trap for every new user-facing builtin, and nothing warns about it. A gate comparing new builtin heads against the linalg convention would catch the next one.`
  - `find_clusters.c still exports exactly one symbol; its machine-precision distance layer (fc_dist_pos, fc_dist_to_point) is private and will be wanted by Classify/Predict NearestNeighbors.`
- Mitigations taken:
  - `A test row pins head == List for both new builtins.`
- Follow-up needed:
  - `DimensionReduce, MultidimensionalScaling, LatentSemanticAnalysis.`

## Activity Log

- `2026-08-14 01:40` Iteration 11: DimensionReduce with all three methods, closing the family's data-in/data-out surface. The three turned out to be ONE algorithm with three ways of forming the matrix to decompose (centred covariance / uncentred Gram / double-centred squared distances), so LSA needed no separate SVD path and MDS reused ml_sym_eigen_desc. Strongest evidence in the family so far: classical MDS on Euclidean distances reproduces PCA exactly to 1e-8, by a completely different route (n x n vs dim x dim), so the agreement is evidence about BOTH rather than self-consistency. Deferred deliberately: a DimensionReducerFunction applicable to new data is a trained model, and that representation belongs with Predict rather than being invented twice.
- `2026-08-14 01:30` Created, and iteration 10 landed: src/ml/pca.{c,h} plus Standardize and PrincipalComponents. Found and fixed a surface bug before it shipped — na_build_matrix returns a visible NDArray, so the first version's `Standardize[d] === {{...}}` was False while Inverse's is True; results now build plain Lists. Canonicalised eigenvector signs so output does not depend on whether LAPACK was linked.

## Reflection

- What went well: `pending`
- What went wrong: `pending`
- Gaps to close: `pending`
- Skill / AGENTS.md updates to propose: `pending`

## Follow-Up

- `Families 3-5: Predict, LearnDistribution (GaussianMixture already reusable via src/ml/gmm.h), Classify.`

## Team Addendum

- `none`
