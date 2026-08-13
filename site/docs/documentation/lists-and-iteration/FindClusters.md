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

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}]
Out[1]= {{1, 2, 3, 1}, {10, 12, 13}, {25}}

In[2]:= FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 4]
Out[2]= {{1, 2, 3, 1}, {10}, {12, 13}, {25}}

In[3]:= FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[4]]
Out[3]= {{1, 2, 3, 1}, {10, 12, 13}, {25}}

In[4]:= FindClusters[{1, 2, 3, 5, 8, 9, 10}, 2]
Out[4]= {{1, 2, 3, 5}, {8, 9, 10}}

In[5]:= FindClusters[{1, a, 3}, 2]
Out[5]= FindClusters[{1, a, 3}, 2]
```

### Options (1)

```mathematica
In[6]:= FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, 3, Method -> "KMeans"]
Out[6]= {{1, 2, 3, 1}, {10, 12, 13}, {25}}
```

## Options & behaviour

### Method

`Method -> m`, or `Method -> {m, subopt -> value}`:

| Method | `Automatic` | `UpTo[n]` | `n` |
|---|:---:|:---:|:---:|
| `"Agglomerate"` (single linkage) | yes | yes | yes |
| `"SpanningTree"` (minimum spanning tree) | yes | yes | yes |
| `"KMeans"` | no | yes | yes |
| `"KMedoids"` | no | yes | yes |
| `"Spectral"` | yes | yes | no |
| `"DBSCAN"` | yes | no | no |
| `"GaussianMixture"` | yes | no | no |
| `"JarvisPatrick"` | yes | no | no |
| `"MeanShift"` | yes | no | no |
| `"NeighborhoodContraction"` | yes | no | no |

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
`ManhattanDistance` and `SquaredEuclideanDistance`. On a line these are monotone
transforms of one another, so **all four give the same partition** for every
method that only ranks distances; the option is accepted for compatibility
rather than offering four behaviours. Any other value leaves the call
unevaluated.

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
  elements or vectors.** The other seven methods read a sorted one-dimensional
  projection and decline rather than run on data for which it is meaningless.
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
- Returns unevaluated for an empty list, a non-`List` argument, a nested list, a
  visible `NDArray`, a non-numeric element, a count that is not a positive
  integer or `UpTo[k]`, or a method incompatible with the count mode.
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

**See also:** [CMYKColor](../../graphics/CMYKColor/), [SquaredEuclideanDistance](../../lists-and-iteration/SquaredEuclideanDistance/), [EditDistance](../../lists-and-iteration/EditDistance/), [Rational](../../arithmetic/Rational/), [Complex](../../arithmetic/Complex/), [List](../../other-advanced/List/), [NDArray](../../linear-algebra/NDArray/), [EuclideanDistance](../../lists-and-iteration/EuclideanDistance/)

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
