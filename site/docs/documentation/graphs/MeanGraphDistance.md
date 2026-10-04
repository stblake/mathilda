# MeanGraphDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MeanGraphDistance[g] gives the mean distance over all ordered pairs of distinct vertices of g: exact for unweighted graphs, a machine real for weighted ones, Infinity unless g is (strongly) connected.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= MeanGraphDistance[PathGraph[4]]
Out[1]= 5/3

In[2]:= MeanGraphDistance[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[2]= Infinity

In[3]:= MeanGraphDistance[Graph[{1},{}]]
Out[3]= MeanGraphDistance[Graph[<1 vertex, 0 edges>]]
```

### Options (1)

```mathematica
In[4]:= MeanGraphDistance[Graph[{1,2,3},{1<->2,2<->3,1<->3}, EdgeWeight->{1,1,5}]]
Out[4]= 1.33333
```

### Applications (2)

Exact for an unweighted graph

```mathematica
In[5]:= MeanGraphDistance[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 5}]]
Out[5]= 2
```

Weighted: a machine real

```mathematica
In[6]:= MeanGraphDistance[Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight -> {1, 4}]]
Out[6]= 3.33333
```

## Implementation notes

**Algorithm.** `builtin_mean_graph_distance` is the mean distance over all ordered pairs of
distinct vertices: `tot / (n(n-1))`, where `tot` is the sum of all pairwise distances taken
from the cached per-source distance summary (`gmet_dist_summary`: multi-source BFS when
unweighted, Dijkstra when weighted). An `O(V+E)` strong-connectivity test short-circuits: a
graph that is not (strongly) connected has mean distance `Infinity`, because some pair is
unreachable.

**Data structures.** The memoised summary's per-vertex distance sums `sum[]`, shared with
`GraphCenter` and `GraphRadius`. Unweighted the result is assembled exactly as an
`Integer`/`Rational` via `make_rational((int64)tot, n(n-1))`; weighted it is a machine `Real`.

**Complexity / limits.** `O(V(V+E))` unweighted, `O(V(E + V log V))` weighted. It is undefined
for a single vertex (`n = 1`) and for the empty graph — both return unevaluated — and
`Infinity` for a non-(strongly-)connected graph; a non-graph argument returns unevaluated.

- *(w)* weight-aware (machine real when weighted).
- Averages over ordered pairs of distinct vertices; exact when unweighted.
- `Infinity` when some pair is unreachable (in particular, unweighted and not
  strongly connected).
- Unevaluated for a single vertex.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gmet_distance.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_distance.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The value is the mean distance over all ordered pairs of distinct vertices. It is exact (an
`Integer` or `Rational`) for an unweighted graph and a machine `Real` for a weighted one.

It is `Infinity` unless the graph is connected (strongly connected, for a directed graph), and
is undefined — left unevaluated — for a graph with a single vertex.
