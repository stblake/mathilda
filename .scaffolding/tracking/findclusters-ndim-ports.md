---
schema_version: 2
title: FindClusters — the eight n-dimensional method ports
slug: findclusters-ndim-ports
status: implementation
source: direct-user-request
owner: Michael Sollami
issue: pending
pull_request: pending
started: 2026-08-13
last_updated: 2026-08-13
blocked_by: none
goal_lock:
  status: active
  stamped: 2026-08-13 23:05
  scope:
    - "src/list/**"
    - "src/ml/**"
    - "src/pack.c"
    - "src/random.c"
    - "src/random.h"
    - "src/core.c"
    - "src/sym_names.c"
    - "makefile"
    - "tests/**"
    - "docs/spec/**"
    - "tools/*.py"
  success_criteria:
    - "AC-1 — All ten methods recover three well-separated blobs in 2-D and 5-D."
    - "AC-2 — For scalar input every method returns exactly what it returned before."
    - "AC-3 — Equal elements are never split across clusters, any method, any dimension."
    - "AC-4 — DistanceFunction changes the partition where Manhattan and Euclidean genuinely differ."
    - "AC-5 — SeedRandom makes the stochastic methods reproducible."
    - "AC-6 — Methods that become quadratic in n-D decline past an explicit cap rather than hanging."
---

# FindClusters — the eight n-dimensional method ports

## Feature Definition

- One-line goal: `Make KMeans, KMedoids, DBSCAN, GaussianMixture, JarvisPatrick, MeanShift, NeighborhoodContraction and Spectral work on vectors, not just scalars.`
- Problem: `find_clusters.c:1790-1793 declines non-scalar input for all but the two MST-based methods. The guard is honest: FcData keeps scalars in double* val consumed through size_t* order and points in a separate double* coord, and the eight non-gap methods reach their data only via fc_sorted_values, which dereferences d->val[d->order[j]] — both NULL for n-D. So these are eight algorithms to write, not a check to relax.`
- Requested by: `direct user request — "work on these", following the NDArray slice in stblake/mathilda#57`

## Current State Study

- Dimension-generic already: `fc_method_gap` (Agglomerate, SpanningTree), which reads only the MST edges.
- 1-D-only helpers the eight share and must stop using: `fc_sorted_values` (`:954`),
  `fc_eps_window` (`:901`), `fc_knn_window` (`:912`), `fc_median_gap` (`:937`),
  `fc_scatter` (`:963`, indexed by *sorted position*), `fc_merge_modes` (`:983`).
- Already in place and free: `d.coord` is populated for every accepted n-D input
  (`:1850-1856`), and the equal-elements fold (`:1914-1934`) is kind-agnostic, so
  `AC-3` holds for the ports without new code.
- The only machine distance is `static fc_sqdist` (`:478`), squared Euclidean only.
- `DistanceFunction` is validated and discarded (`:1658`); `FcOpts` has no field for it.

## Scope

### In Scope

- A generic distance entry point, and making `DistanceFunction` actually apply
- An n-D length scale replacing the 1-D `fc_median_gap` defaults
- Porting all eight methods, and relaxing the guard as each becomes correct
- Explicit caps for the methods that become quadratic
- Tests, gates, docs, changelog

### Out Of Scope

- `Classify`, `Predict`, `LearnDistribution`, `DimensionReduce` — separate families
- Approximate neighbour structures (k-d / ball trees); brute force matches Prim's cost
- `CriterionFunction` — still accepted and ignored; making it real is its own piece

## Success Criteria

- `AC-1` — all ten methods recover three well-separated blobs in 2-D and 5-D
- `AC-2` — scalar input unchanged for every method
- `AC-3` — equal elements never split
- `AC-4` — `DistanceFunction` changes the partition where the metrics differ
- `AC-5` — `SeedRandom` reproduces the stochastic methods
- `AC-6` — quadratic methods decline past a cap rather than hanging

## Tests

- Must-have: blobs in 2-D and 5-D per method; scalar regression per method;
  `{{1},{2},{100}}` (POINT with dim 1); duplicated points; DistanceFunction
  difference; SeedRandom reproducibility; cap decline.
- Framework: the CMake suite. A new file needs **both** `add_executable` and
  `add_test` in `tests/CMakeLists.txt`; nothing is auto-discovered.

## Task List

- [ ] `prerequisites: fc_dist, n-D length scale, assign[] contract, DistanceFunction` | `independent` | `in-progress`
- [ ] `MeanShift + NeighborhoodContraction` | `depends on: prerequisites` | `pending`
- [ ] `KMeans` | `depends on: prerequisites` | `pending`
- [ ] `DBSCAN` | `depends on: prerequisites` | `pending`
- [ ] `JarvisPatrick` | `depends on: prerequisites` | `pending`
- [ ] `KMedoids` | `depends on: prerequisites` | `pending`
- [ ] `Spectral` | `depends on: KMeans` | `pending`
- [ ] `GaussianMixture` | `depends on: prerequisites` | `pending`
- [ ] `tests, gates, docs, changelog` | `depends on: all` | `pending`

## Checkpoints

- [x] start | completed: `2026-08-13 23:05`
- [x] spec / plan created | completed: `2026-08-13 23:05`
- [x] implementation started | completed: `2026-08-13 23:05`
- [ ] implementation complete | completed: `pending`
- [ ] critic pass | completed: `pending`
- [ ] risk-register reviewed | completed: `pending`
- [ ] feature validated | completed: `pending`
- [ ] PR created | completed: `pending`
- [ ] closeout complete | completed: `pending`

## Decisions

- `Ported methods write assign[] directly rather than going through fc_scatter, which is indexed by sorted position and has no n-D meaning. fc_method_gap already does this, so the ports converge on its contract rather than inventing a third one.`
- `The n-D length scale is the median MST edge weight, already computed in d->gap, so it costs nothing and is the natural analogue of the 1-D median adjacent gap.`

## Risks And Unknowns

- `The shared refactor touches code all ten methods use, including the 1-D paths that work today. AC-2 must be written and passing BEFORE the refactor, not after, or a regression is invisible.`
- `Spectral's cut rule has no n-D analogue, so it depends on the KMeans port — more coupling than the task list suggests.`
- `KMeans' empty-cluster repair depends on a "cluster is a contiguous sorted run" invariant that does not exist in n-D.`

## Risk Register

| ID | Category | Likelihood | Impact | Mitigation | Residual | Owner |
|---|---|---|---|---|---|---|
| RISK-1 | Correctness regression | med | high | The shared refactor touches the working 1-D paths. `AC-2` pins every method's scalar output before any refactor lands. | Low if the pins go first. | Michael Sollami |
| RISK-2 | Denial of service (self-inflicted) | med | med | Methods that become O(n²) with no cap could hang the kernel, and there is no cooperative abort — the only exit is killing the process. Explicit caps, declining rather than hanging. | Accepted: a decline is visible and safe. | Michael Sollami |

## Test Results

- Command: `pending`
- Outcome: `pending`

## Proofs Of Completion

- Completion timestamp: `pending`

## Tech Debt Review

- Potential tech debt introduced:
  - `pending`
- Existing tech debt noticed:
  - `CriterionFunction accepted and ignored (find_clusters.c:1677)`
  - `CosineDistance and HammingDistance implemented in distance.c but not in FindClusters' accepted names`
- Mitigations taken:
  - `pending`
- Follow-up needed:
  - `pending`

## Activity Log

- `2026-08-13 23:05` Created. Follows the NDArray slice (stblake/mathilda#57); the user asked to work through the ML families, starting by completing FindClusters.

## Reflection

- What went well: `pending`
- What went wrong: `pending`
- Gaps to close: `pending`
- Skill / AGENTS.md updates to propose: `pending`

## Follow-Up

- `The other four Wolfram ML families, cheapest first: DimensionReduce + PCA, then Predict + LinearRegression, then LearnDistribution, then Classify. See the mathilda-ml-followups memory note.`

## Team Addendum

- `none`
