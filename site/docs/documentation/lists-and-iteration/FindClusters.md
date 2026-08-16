# FindClusters

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindClusters[list]`**

Partitions a 1D numeric list into clusters of nearby elements, as a list of lists. Clusters appear in order of the first occurrence of a member; elements keep their input order.

**`FindClusters[list, n]`**

Gives exactly n clusters, capped at the number of distinct values.

**`FindClusters[list, UpTo[n]]`**

Gives at most n clusters, and fewer when the data suggests fewer.

**`FindClusters[list, spec, Method -> m]`**

Uses algorithm m: Agglomerate, SpanningTree, KMeans, KMedoids, Spectral, DBSCAN, GaussianMixture, JarvisPatrick, MeanShift or NeighborhoodContraction. KMeans and KMedoids require a count; the density methods require Automatic. Options: Method, DistanceFunction (Automatic, EuclideanDistance, ManhattanDistance or SquaredEuclideanDistance -- all equivalent in 1D), CriterionFunction and PerformanceGoal (accepted, no effect). Returns unevaluated for a non-numeric element, an empty list, a method incompatible with the count mode, or a list too large for the chosen method (Spectral above 2000 elements, MeanShift and NeighborhoodContraction above 4000, both being quadratic).

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= FindClusters[{{0, 0}, {0, 11}, {8, 6}}, 2]
Out[1]= {{{0, 0}}, {{0, 11}, {8, 6}}}
```

### Scope (5)

```mathematica
In[2]:= FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}]
Out[2]= {{1, 2, 3, 1}, {10, 12, 13}, {25}}

In[3]:= FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 4]
Out[3]= {{1, 2, 3, 1}, {10}, {12, 13}, {25}}

In[4]:= FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[4]]
Out[4]= {{1, 2, 3, 1}, {10, 12, 13}, {25}}

In[5]:= FindClusters[{1, 2, 3, 5, 8, 9, 10}, 2]
Out[5]= {{1, 2, 3, 5}, {8, 9, 10}}

In[6]:= FindClusters[{1, a, 3}, 2]
Out[6]= FindClusters[{1, a, 3}, 2]
```

### Options (2)

```mathematica
In[7]:= FindClusters[{{0, 0}, {0, 11}, {8, 6}}, 2, DistanceFunction -> ManhattanDistance]
Out[7]= {{{0, 0}, {0, 11}}, {{8, 6}}}

In[8]:= FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> "KMeans"]
Out[8]= {{1, 2, 3, 1}, {10, 12, 13}, {25}}
```

## Options & behaviour

### Method

`Method -> m`, or `Method -> {m, subopt -> value}`:

The last column is the dimensionality: `vectors` means the method clusters
n-dimensional points, `scalars` means it still reads the sorted one-dimensional
projection and declines vector input.

| Method | `Automatic` | `UpTo[n]` | `n` | data |
|---|:---:|:---:|:---:|:---:|
| `"Agglomerate"` (single linkage) | yes | yes | yes | vectors |
| `"SpanningTree"` (minimum spanning tree) | yes | yes | yes | vectors |
| `"KMeans"` | no | yes | yes | vectors |
| `"KMedoids"` | no | yes | yes | vectors |
| `"Spectral"` | yes | yes | no | vectors |
| `"DBSCAN"` | yes | no | no | vectors |
| `"GaussianMixture"` | yes | no | no | vectors |
| `"JarvisPatrick"` | yes | no | no | vectors |
| `"MeanShift"` | yes | no | no | vectors |
| `"NeighborhoodContraction"` | yes | no | no | vectors |

- `Automatic` (the default) currently means `"Agglomerate"`. It is not a
  criterion-driven search across methods.
- **`"Agglomerate"` and `"SpanningTree"` are the same computation in 1D**, and
  share one implementation. Single-linkage clustering equals cutting the widest
  edges of the minimum spanning tree, and on a line the spanning tree is the
  sorted adjacency chain.
- With `Automatic`, the gap methods cut every sorted-adjacent gap wider than
  three times the median gap; uniform spacing therefore gives a single cluster.
- Suboptions: `"NeighborhoodRadius"` (`"DBSCAN"`, `"MeanShift"`,
  `"NeighborhoodContraction"`), `"MinPoints"` (`"DBSCAN"`), `"NeighborCount"`
  (`"JarvisPatrick"`). An unrecognised suboption leaves the call unevaluated.
- `"DBSCAN"` returns noise points as singleton clusters, so the result is always
  a partition of the input.
- `"KMedoids"` is Voronoi iteration, not PAM, so it can settle on a worse
  partition than a swap-based search would.
- `"Spectral"` builds an n×n similarity matrix and declines above 2000 elements.
  In 1D its embedding largely recovers the sorted order, so it rarely differs
  much from `"Agglomerate"`.

### DistanceFunction

accepts `Automatic`, `EuclideanDistance`,
`ManhattanDistance` and `SquaredEuclideanDistance`, and **selects the metric the
spanning tree is weighted by**. Takes a metric, not a string: `EuclideanDistance`,
not `"EuclideanDistance"` — the reverse of `Method`, which takes strings. Any
other value leaves the call unevaluated, so a typo cannot silently cluster by the
wrong metric.

Two things follow from the arithmetic rather than from choice:

- **On a line all four agree.** Every accepted metric is a monotone transform of
  `|a - b|` there, so they induce the same ordering of gaps and the same
  partition. The option is therefore accepted and has no effect on scalar input,
  whose weights are plain exact differences.
- **`EuclideanDistance` and `SquaredEuclideanDistance` always agree**, in any
  dimension. Squaring is monotone on non-negatives so it preserves edge ranking,
  and the `Automatic` threshold compares against a multiple of the median, where
  `d > 3 median(d)` exactly when `d² > 9 median(d²)`. They are one implementation
  for that reason — which also avoids taking a square root of an exact value,
  since `Sqrt[2]` is irrational and an exact ordering that must compare
  irrationals is a far harder problem than clustering needs.

`ManhattanDistance` is the one that genuinely differs above one dimension:

Squared Euclidean makes those edges 121, 100 and 89, so the tree is `{89, 100}`
and cutting its heaviest edge isolates `{0, 0}`. Manhattan makes them 11, 14 and
13, so the tree is `{11, 13}` and cutting its heaviest edge isolates `{8, 6}`
instead.

**CriterionFunction** and **PerformanceGoal** are accepted and currently have no
effect.

## Algorithm

FindClusters[list] / [list, n] / [list, UpTo[n]] -- partition a 1D numeric list into clusters of nearby elements.

```text
  FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}]  -> {{1, 2, 3, 1}, {10, 12, 13}, {25}}
  FindClusters[{1, 2, 3, 5, 8, 9, 10}, 2]     -> {{1, 2, 3, 5}, {8, 9, 10}}
```

WE DO NOT REPRODUCE MATHEMATICA'S OUTPUT, and this is not a shortfall to be fixed later -- it is not reachable. Mathematica auto-selects a distance function and preprocesses the data by rules it does not publish; its own messages name the choice ("... method Agglomerate and distance function HammingDistance ..."). The visible consequence is that even single-linkage with an explicit count disagrees with the textbook algorithm:

```text
  FindClusters[{1,4,9,16,25,36}, 3, Method->"Agglomerate"] -> {{36},{1,4,9},{16,25}}
```

on gaps 3,5,7,9,11 it cut 11 and 7 rather than the two largest, 11 and 9. Its Automatic cluster count is likewise not a gap rule -- {1,4,9,16,25,36} (linear gap growth) splits three ways while {1,2,4,8,16,32,64} (geometric, far more extreme gaps) gives one. So: this file implements the textbook algorithm for each named method with the semantics stated below, and the acceptance table in tests/test_list.c is the specification.

### Semantics

```text
  Distance          |a - b| on the line. Never via Abs -- see EXACTNESS.
  Cluster order     by first occurrence of any member in the input.
                    Mathematica's is unstable (Automatic and n=3 disagree on
                    the same input); first occurrence matches Gather.
  Element order     input order, within every cluster.
  Count             three modes, see FcCountMode. `n` is capped at the
                    distinct-value count, never exceeded.
```

### Exactness

The sorted order and the fixed-count gap selection are computed on the ELEMENTS THEMSELVES via list_numeric_cmp, so exact input is ordered exactly: FindClusters[{1/10^25, 2/10^25, 1}, 2] works. Nearest cannot do this, because it routes its distance through Abs, and builtin_abs declines on a rational with a bigint component (complex.c:418-421). A sorted 1D pass has no need of Abs at all -- the gap between sorted neighbours is non-negative by construction -- so the landmine is avoided rather than inherited. Rewriting a gap as Abs[b - a] would reintroduce it, and an acceptance row fails if anyone does.

The Automatic threshold and the inherently-inexact methods (KMeans, the density family, GaussianMixture, Spectral) work on a machine-double projection instead, which is correct: a mean, a kernel and an eigenvector are inexact by definition, and the Automatic threshold carries a fitted constant.

### Cost

Every method is O(n log n), dominated by the initial sort, except Spectral, which is O(n^2) in memory and declines above FC_SPECTRAL_MAX_N. The density family stays linear only because the array is sorted -- the same neighbourhood queries in general dimension are the expensive part of those algorithms.

### Scope

1D numeric lists. Mathematica also clusters symbolic elements as nominal features (FindClusters[{1, a, 3}, 2] gives {{a}, {1, 3}}); we decline. The rule forms, Association input, CriterionFunction, PerformanceGoal, Weights and the FeatureX options are not implemented.

## Performance

Measured on arm64 Darwin at commit `2dea9cc05`.

| case | n | time |
|---|---:|---:|
| Agglomerate (default) | 1,000 | 545 us |
| Agglomerate (default) | 10,000 | 6.5 ms |
| Agglomerate (default) | 100,000 | 96.9 ms |
| KMeans, 5 clusters | 1,000 | 423 us |
| KMeans, 5 clusters | 10,000 | 4.5 ms |
| KMeans, 5 clusters | 100,000 | 55.9 ms |
| 2-D points | 500 | 101.0 ms |
| 2-D points | 1,000 | 397.8 ms |
| 2-D points | 2,000 | 1.55 s |

## Implementation notes

- `Protected`.
- The result is a list of lists. Clusters appear in order of the first
  occurrence of any member; elements keep their input order within a cluster.
- **Element kinds and the metric each uses.** All elements must be of one kind;
  a mixture leaves the call unevaluated.
  - Real numbers -- absolute difference on the line. This is the original path
    and is unchanged.
  - Equal-length numeric vectors (`{{2.5, 3.1}, {5.9, 3.4}, ...}`) -- squared
    Euclidean distance. Ranking on the square rather than the root is what keeps
    the partition exact for exact input, since squaring is monotone on
    non-negatives and `SquaredEuclideanDistance` is rational for rational input.
  - Colours, whose arguments are coordinates: `RGBColor[r, g, b]` is a point in
    its own space and clusters like a 3-vector. Every element must share one
    colour head; `RGBColor` and `GrayLevel` together are declined.
  - Strings -- `EditDistance`, the Levenshtein distance.
  - `Rational` and `Complex` are *numbers*, not points: `{1/2, 1/3, 10}` clusters
    as three reals rather than as pairs of coordinates, and complex input is
    declined as before.
- **In one dimension the spanning tree IS the sorted adjacency chain**, so the
  general implementation reduces exactly to the original algorithm there and
  one-dimensional results are unchanged. Above one dimension the tree is a real
  minimum spanning tree built by Prim's algorithm over exact distances.
- **Only `Automatic`, `"Agglomerate"` and `"SpanningTree"` accept non-numeric
  elements.** For *vectors* the list is longer — see "Which methods accept vectors"
  below. The methods that still read a sorted one-dimensional projection decline
  rather than run on data for which it is meaningless.
- **Two spanning-tree builders, and therefore two ceilings.** Points whose
  every coordinate is already a machine number take a double-precision Prim,
  since such input carries no precision that exact arithmetic would preserve;
  points with a `Rational`, bigint or MPFR coordinate keep the exact builder,
  which is the only one that can order them correctly. Distinctness is decided by
  comparing the *elements* exactly on both paths, so `{{1/3, 1/7}, {1/3, 1/7 +
  1/10^20}}` stays together either way.
  - machine points: **capped at 20000**. Measured 6.8 ms at 2000, 0.17 s at
    10000, 0.81 s at 20000.
  - exact points and strings: **capped at 2000**, at roughly 1.5 s there.
  - Both are quadratic -- a sort makes neighbourhood work linear only on a line
    -- so a larger input is declined rather than silently taking minutes. One
    dimension has no cap and does 10^6 elements in about 2.3 s.
- `n` is capped at the number of distinct values — no method can separate two
  equal elements, so `FindClusters[{1, 2, 3}, 5]` gives three clusters and
  `FindClusters[{7, 7, 7, 7}, 3]` gives one.
- The distinction between the two count forms is real:
  `FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 4]` forces a fourth cluster by
  splitting `{10, 12, 13}`, while `UpTo[4]` returns the natural three.
- **Exact input is ordered exactly.** The sort and the fixed-count gap selection
  compare the elements themselves, so rationals with bigint components work:
  `FindClusters[{1/10^25, 2/10^25, 1}, 2]` is
  `{{1/10000000000000000000000000, 1/5000000000000000000000000}, {1}}`.
- **A visible `NDArray` is accepted**, at rank 1 (scalars) or rank 2 (points), and
  gives the same partition as the equivalent `List`:
  `FindClusters[NDArray[{1., 2., 10.}]]` is `{{1., 2.}, {10.}}`. It is
  materialised on the way in rather than read as a buffer, because the exact
  spanning tree, the exact boundary set and the emitted clusters all work on the
  input elements themselves — and that exactness is why one-dimensional answers
  are exact. Nothing is given up in speed: every value in an `NDArray` is machine
  by construction, so the machine spanning-tree builder runs, the same one a
  machine `List` takes. A *packed* list needs no special handling — the
  transparency gate materialises it before the builtin sees it, which is why only
  the visible surface needed a guard. Rank 3 and above still decline: a
  list-valued component is not a coordinate.
- **Which methods accept vectors.** `Agglomerate`, `SpanningTree`, `MeanShift`,
  `NeighborhoodContraction`, `KMeans`, `DBSCAN`, `JarvisPatrick`, `KMedoids`,
  `Spectral` and `GaussianMixture` cluster n-dimensional points — **all ten**. The
  only input still declined above one dimension is a list of strings, which has no
  coordinates; the gap methods handle those through edit distances in the tree and return unevaluated for vector input.
  That is not a conservative guard — those five reach their data through the sorted
  projection, which does not exist off a line — and the list grows as each is
  ported. The ported methods additionally decline string input, having no
  coordinates to move; the gap methods still handle strings through edit distances
  in the tree.
- **`KMeans` above one dimension is a second implementation, not the same one.**
  The one-dimensional kernel seeds at quantiles of the sorted distinct values; the
  n-dimensional one seeds *farthest-first* — the point nearest the centroid, then
  repeatedly the point farthest from everything already chosen. Both are
  deterministic and need no random seeding. Two consequences worth relying on:
  - the result **does not depend on the order the points were given in**, because
    the first centre is chosen from a property of the set rather than from the
    first element; and
  - every centre starts as a distinct data point, so empty clusters are rare. One
    that does empty is reseeded at the point currently farthest from its own
    centre, which is the order-free form of the one-dimensional repair.

  The two paths agree on the same data written both ways
  (`FindClusters[{1, 2, 3, 10, 11, 12, 25}, 3, Method -> "KMeans"]` and the same
  numbers as `{{1}, {2}, …}`), and unifying them would move the one-dimensional
  answers, so it stays a deliberate change rather than a tidy-up.
- **`"JarvisPatrick"` needs enough points per group for its `NeighborCount`.** The
  default is 5, and a 5-nearest-neighbour list cannot fit inside a 4-point group —
  it must reach into a neighbouring one, which then links them. So
  `FindClusters[fourPointBlobs, Method -> "JarvisPatrick"]` returns *one* cluster
  while `"NeighborCount" -> 3` returns three. That is the algorithm on data too
  small for its default rather than a defect, and the one-dimensional kernel does
  the same (a `k` clamped to `n - 1` puts every point in one window).
- **`"JarvisPatrick"` above one dimension is a second implementation**, like
  `KMeans` and unlike `DBSCAN`. The 1-D kernel counts shared neighbours as the
  overlap of two *contiguous windows* and links only adjacent sorted pairs, where
  the n-D one takes a true set intersection over all pairs — the textbook rule.
  The two differ on at least one checked case
  (`{1, 2, 10, 12, 3, 1, 13, 25}` at `"NeighborCount" -> 2`), so adopting the
  general rule on a line would move an existing answer and is left as its own
  deliberate change.
- **`"KMedoids"` carries a tighter ceiling than `"KMeans"`,** because the two differ
  by a complexity class rather than a constant. A mean is `O(n·dim)`; a medoid
  search compares every member against every other member of its cluster and is
  `O(n²·dim)` however small `k` is. So the same 2000×10 input that `"KMeans"`
  clusters is *declined* by `"KMedoids"`.
- **`"KMedoids"` above one dimension is a second implementation** (quantile seeding
  on a line, farthest-first off it), and the two **disagree** on
  `{1, 2, 3, 10, 11, 12, 25}` at `k = 3`: the 1-D kernel gives
  `{{1,2,3}, {10}, {11,12,25}}` and the n-D kernel `{{1,2,3}, {10,11,12}, {25}}`.
  The n-D answer is strictly better by the method's own objective — total distance
  from each member to its cluster's medoid — **4 against 16**. Both are documented
  as they behave; adopting farthest-first on a line would improve the 1-D answer and
  move an existing one, so it is left as its own deliberate change.
- **`"Spectral"` above one dimension reads the graph's connected components before
  its Fiedler vector, and that is not an optimisation.** The multiplicity of the
  zero eigenvalue of the normalised Laplacian equals the number of connected
  components, so on a well-separated sample the null space is multi-dimensional and
  a single "leading non-trivial eigenvector" is defined only up to a rotation
  inside it — power iteration then lands on whichever member the seed favours,
  separating two groups at best. This is the *easy* case, not a corner case: better
  separation means smaller cross-affinity, and `exp(-(60/1.5)²/2)` is zero in double
  precision, not merely small. Components are therefore read off with union-find.
  Asked (via `UpTo[k]`) for fewer clusters than there are components, the nearest
  components merge first, ranked by the spanning tree — the spectrum itself cannot
  choose, since every partition respecting components has zero cut.
- **`"Spectral"` above one dimension is a second implementation, and the threshold
  is why.** The n-D kernel counts embedding jumps wider than three times the **mean**
  jump; the 1-D kernel counts data gaps wider than three times the **median**. Each
  is right in its own domain — a good embedding collapses within-cluster distance,
  driving the median to ~0 so a median threshold accepts nearly every jump, while on
  raw data gaps it is the mean that misleads. Routing scalars through the n-D kernel
  under-counts (`{{1,2,3},{10,11,12},{25}}` becomes `{{1,2,3,10,11,12},{25}}`), so
  the two paths stay separate.
- **`"GaussianMixture"` fits a full covariance, so it needs more points per
  component than dimensions.** A component costs `dim` means plus `dim(dim+1)/2`
  covariance entries, and `k_max` is capped at `n/(dim + 1)` — the minimum for a
  non-degenerate scatter matrix. So nine points in five dimensions correctly return
  *one* cluster (`k_max` falls to 1), while twenty-four points per blob recover all
  three. This is BIC and the parameter count doing their job, not a defect.
- **A singular component is modelled, not refused.** Identical points give a zero
  scatter matrix and collinear points a covariance with a zero eigenvalue. The
  variance floor is added as a **ridge** on the covariance diagonal — the matrix
  reading of the scalar "clamp the variance up to a minimum", chosen because clamping
  every direction properly means an eigendecomposition per component per iteration
  while a ridge is one add and can only overshoot the bound. A component whose full
  covariance still will not factorise falls back to diagonal rather than failing the
  fit.
- **`"GaussianMixture"`'s kernel is the one clustering algorithm that does not live
  in `find_clusters.c`.** The EM fit is in `src/ml/gmm.c` behind a buffer-level API
  (a row-major `n × dim` array, not a `FindClusters` internal), because it is the
  first kernel with two real consumers — this `Method` and the `LearnDistribution`
  learner of the same name. `find_clusters.c` exports exactly one symbol, so a second
  builtin could not otherwise reach it.
- **`KMeans` and `KMedoids` require a count.** Neither accepts the `Automatic`
  form, in any dimension, so `FindClusters[data, Method -> "KMeans"]` returns
  unevaluated; `UpTo[k]` is the data-driven form they do take, and a bare `k` is
  obeyed exactly.
- **The n-dimensional `KMeans` ceiling is on `n × k × dim`, not on `n`.** Lloyd's
  algorithm is linear in `n` — cheaper than the spanning tree already built for the
  same input — so an `n` cap would refuse work just paid for: 20000 machine points
  in ten dimensions at `k = 3` takes **0.60 s**, and at `UpTo[8]` the same. What is
  declined is a `k` on the order of thousands at that size, where the assignment
  scan becomes a hundred near-quadratic passes over the whole sample.
- Returns unevaluated for an empty list, a non-`List` argument, a nested list, an
  `NDArray` of rank 3 or more, a non-numeric element, a count that is not a
  positive integer or `UpTo[k]`, or a method incompatible with the count mode.
- **Symbolic elements are declined**, not clustered. Mathematica treats them as
  nominal features and answers `FindClusters[{1, a, 3}, 2]` with `{{a}, {1, 3}}`;
  this implementation is numeric-only.
- **The output is not expected to match Mathematica's.** Mathematica
  auto-selects a distance function and preprocesses the data by unpublished
  rules, so even single-linkage with an explicit count differs, and its
  `Automatic` cluster count is an internal index rather than a gap rule. The
  algorithms below are the textbook ones with the semantics stated here.

**Attributes:** `Protected`.

## References

**See also:** [CMYKColor](../../graphics/CMYKColor/), [SquaredEuclideanDistance](../../lists-and-iteration/SquaredEuclideanDistance/), [EditDistance](../../lists-and-iteration/EditDistance/), [Rational](../../arithmetic/Rational/), [Complex](../../arithmetic/Complex/), [NDArray](../../linear-algebra/NDArray/), [List](../../other-advanced/List/), [LearnDistribution](../../machine-learning/LearnDistribution/)

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_findclusters_distance.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findclusters_distance.c)
- Tests: [`tests/test_findclusters_ndim.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findclusters_ndim.c)
- Tests: [`tests/test_findclusters_scalar_pin.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findclusters_scalar_pin.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
