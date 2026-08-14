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

- [x] `DistanceFunction wired through FcOpts and both MST builders` | `independent` | `done (iteration 2)`
- [x] `prerequisites: fc_points, fc_dist_to_point, fc_dist_pos, fc_scale_ndim, union-find merge` | `independent` | `done (iteration 3, landed with their first consumer)`
- [x] `MeanShift + NeighborhoodContraction` | `depends on: prerequisites` | `done (iteration 3)`
- [x] `KMeans` | `depends on: prerequisites` | `done (iteration 4)`
- [x] `DBSCAN` | `depends on: prerequisites` | `done (iteration 5) -- REPLACED the 1-D kernel; pins proved the general rule answer-preserving`
- [x] `JarvisPatrick` | `depends on: prerequisites` | `done (iteration 6) -- two paths; general rule moves a list_tests answer`
- [x] `KMedoids` | `depends on: prerequisites` | `done (iteration 7) -- reuses fc_lloyd_ndim behind a flag; second tighter cap; 1-D/n-D disagree and n-D is better (cost 4 vs 16)`
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
  - `NeighborhoodContraction cannot separate an outlier from a tie-heavy 4-point set: a zero median edge weight makes the scale fall back to range/(n-1), and the flat kernel default radius of 3*scale is then EXACTLY the range, putting the outlier precisely on the inclusion boundary. Pre-existing 1-D behaviour, identical in both dimensionalities. Fixing it moves a pinned answer, so it needs its own deliberate change.`
  - `JarvisPatrick's 1-D rule is weaker than the textbook one: shared neighbours as a contiguous-window overlap, linkage restricted to adjacent sorted pairs. The n-D kernel does the real thing. Adopting it on a line changes FindClusters[{1,2,10,12,3,1,13,25}, Method -> {"JarvisPatrick", "NeighborCount" -> 2}] from {{1,2,1},{10,12,13},{3},{25}} to {{1,2,3,1},{10,12,13},{25}}. Probably an improvement; needs its own deliberate change.`
  - `1-D KMedoids sits in a worse local optimum than the n-D kernel on the same data: FindClusters[{1,2,3,10,11,12,25}, 3, Method -> "KMedoids"] gives cost 16 where farthest-first gives cost 4. Adopting farthest-first on a line would improve it and move a pinned answer -- its own deliberate change.`
  - `EVERY clustering kernel is static in find_clusters.c: 55 static functions, ONE exported symbol (builtin_find_clusters). ClusteringComponents, ClusteringTree, Dendrogram, NearestNeighborGraph do not exist and cannot be built without extraction. GaussianMixture is needed by BOTH FindClusters and LearnDistribution, so it is the forcing function. distance.c exports an EXACT metric layer (distance_squared_euclidean etc.) but the MACHINE-precision layer (fc_dist_pos, fc_dist_to_point) is private -- that is the layer every ML family needs for speed, and the duplication risk for families 2-5.`
  - `CosineDistance and HammingDistance implemented in distance.c but not in FindClusters' accepted names`
  - `KMeans now has two Lloyd implementations, one per dimensionality. Unifying them onto farthest-first would move the pinned 1-D answers, so it is a deliberate behaviour change and not a refactor. The n-D suite asserts the two agree on the same data written both ways, which is what would catch an accidental merge.`
  - `The KMeans work budget refuses k == n, which actually converges in a single iteration since farthest-first hands every point its own centre. Conservative on purpose -- modelling the iteration count per k is more machinery than the degenerate case is worth -- but it is a false refusal.`
- Mitigations taken:
  - `pending`
- Follow-up needed:
  - `pending`

## Activity Log

- `2026-08-14 00:50` Iteration 7: KMedoids ported by reusing fc_lloyd_ndim behind a medoid flag rather than duplicating 150 lines -- farthest-first already seeds with real data points, which is what a medoid needs. Added a SECOND, tighter ceiling because the medoid search is a different complexity class (O(n^2 dim) vs O(n dim)); the same 2000x10 input now separates the two methods, which is the cleanest proof the cap is real. The two dimensionalities DISAGREE here and the n-D answer is strictly better by the method's own objective, 4 vs 16 -- so the 1-D kernel sits in a 4x worse local optimum, recorded as tech debt rather than silently fixed. Negative result also recorded: 26 constructed cases failed to make KMedoids and KMeans produce different PARTITIONS, so no test claims they do.
- `2026-08-14 01:25` Iteration 6: JarvisPatrick ported with real k-NN lists and set intersection. Kept TWO paths, and finding out why corrected the procedure from iteration 5: routing scalars through the general kernel leaves all 22 pins passing -- the DBSCAN criterion -- yet moves a list_tests answer at NeighborCount 2, because the 1-D kernel counts shared neighbours as the overlap of two contiguous windows and links only adjacent sorted pairs. The unify test is therefore 'does any COVERED 1-D answer move', not 'do the pins pass'; the pin file is a floor, not the whole contract. Also learned that this method's default NeighborCount 5 cannot fit inside the 4-point blobs the rest of the suite uses, so its blob test needs 8-per-blob -- asserting the merge would have pinned a data-sizing artefact. Measured DBSCAN's ceiling too, which the previous entry had asserted without a number: 1.73 s at 20000x10.
- `2026-08-14 00:55` Iteration 5: DBSCAN ported, and unified onto ONE kernel rather than two. The 1-D rule linked only adjacent sorted pairs where DBSCAN specifies any eps-close pair through a core point; at the default MinPoints 2 those provably coincide, above it the argument runs out, so I put the question to the pin suite instead of reasoning about it -- all 22 scalar pins pass through the general kernel, so the 1-D body and fc_eps_window were deleted. This is the pin suite paying for itself: it is what let DBSCAN unify and KMeans stay split, on evidence rather than taste.
- `2026-08-14 00:05` Iteration 4: KMeans ported. NOT unified with the 1-D kernel, and the contrast with iteration 3 is the point: MeanShift unified because the n-D form provably reproduced the 1-D answers, whereas quantile and farthest-first initialisation settle in different local optima, so routing one through the other would move pinned answers. Two paths, with a test that asserts they agree on the same data written both ways. Farthest-first chosen for three properties, not by default: deterministic (no RandomVariate needed yet), independent of input order (the first centre comes from the centroid, a property of the SET, not from index 0), and every centre a distinct data point so empty clusters are rare. Corrected my own first cap: I copied the shift methods' n-cap by analogy, which was the wrong SHAPE — Lloyd is linear in n and cheaper than the spanning tree already built, so an n-cap declines work just paid for. Rewrote it as a budget on n*k*dim; 20000x10 at k=3 now runs in 0.60 s where the first version refused it.
- `2026-08-13 23:35` Iteration 3: MeanShift and NeighborhoodContraction ported. Unified rather than duplicated, on three observations: d->val is already the dim-1 layout of d->coord; the median MST edge weight generalises the median adjacent gap exactly, and the mean generalises the zero-median fallback, because on a line the MST is the sorted chain whose gaps sum to the range; and union-find merging reproduces the adjacent-difference pass on a line. fc_merge_modes deleted as irreducibly 1-D. Found and recorded a pre-existing NeighborhoodContraction quirk on tie-heavy 4-point input.
- `2026-08-13 23:20` Iteration 2: DistanceFunction made real. It was validated and discarded, which is harmless on a line (all four metrics are monotone transforms of |a-b| there) and a wrong answer above one dimension. Euclidean and SquaredEuclidean collapsed to one case on purpose — monotone, same threshold test, and it avoids an exact square root. The threshold factor now tracks the metric rather than the kind, since Manhattan weights are linear and the squared factor would have been nine times too large.
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
