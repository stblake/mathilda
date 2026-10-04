# GraphPeriphery

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphPeriphery[g] gives the vertices of g with maximum eccentricity; {} unless g is (strongly) connected.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= GraphCenter[PathGraph[5]]
Out[1]= {3}

In[2]:= GraphPeriphery[PathGraph[5]]
Out[2]= {1, 5}

In[3]:= GraphCenter[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= {}

In[4]:= GraphPeriphery[Graph[{1<->2,3<->4}]]
Out[4]= {}
```

### Options (1)

```mathematica
In[5]:= GraphCenter[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}]]
Out[5]= {1}
```

### Applications (3)

The two endpoints are farthest out

```mathematica
In[6]:= GraphPeriphery[PathGraph[{1, 2, 3, 4, 5}]]
Out[6]= {1, 5}
```

A vertex-transitive graph: every vertex is peripheral

```mathematica
In[7]:= GraphPeriphery[CycleGraph[5]]
Out[7]= {1, 2, 3, 4, 5}
```

All eccentricities equal 1

```mathematica
In[8]:= GraphPeriphery[CompleteGraph[4]]
Out[8]= {1, 2, 3, 4}
```

## Implementation notes

**Algorithm.** `builtin_graph_periphery` gives the vertices of maximum
eccentricity — the "rim" of the graph — where a vertex's eccentricity is the
greatest distance from it to any other vertex. It shares the `extremal`
dispatcher with `GraphDiameter`/`GraphRadius`/`GraphCenter`/`MeanGraphDistance`:
an all-pairs distance summary (`gmet_dist_summary`) yields each vertex's
eccentricity, and the periphery is the set achieving the maximum. For a weighted
graph a vertex that cannot reach every other has eccentricity `Infinity` (as in
the Wolfram Language), and ties are compared with a small relative tolerance so
two equal weight-sums taken in different orders match. If `g` is not (strongly)
connected the periphery is `{}`.

**Data structures.** A `GmetCSR` adjacency plus the cached `GmetDistSummary`
(per-vertex eccentricity `ecc[]`, reach counts, distance sums). The matching
vertices are collected through `gmet_vertex_subset`, and the answer is cached on
the graph node (`gmet_cache_put`).

**Complexity / limits.** Dominated by the all-pairs distance computation (BFS per
source unweighted, Dijkstra per source weighted). Integer distances for an
unweighted graph. `{}` unless the graph is connected.

- *(w)* weight-aware.
- Unweighted and not strongly connected (not connected, if undirected): `{}`.
- Weighted: taken over the weighted eccentricities, so a weighted digraph that
  is not strongly connected can have a non-empty center
  (*reverse-engineered*).
- Reduced from the cached MS-BFS per-source summary (see `GraphDistanceMatrix`).

**Attributes:** `Protected`.

## References

**See also:** [GraphCenter](../../graphs/GraphCenter/), [GraphDistanceMatrix](../../graphs/GraphDistanceMatrix/)

- Source: [`src/graph/gmet_distance.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_distance.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The periphery is the set of vertices whose eccentricity equals the graph
diameter — the vertices "on the rim". Dually, `GraphCenter` picks those of
minimum eccentricity (the radius).

On a vertex-transitive graph (a cycle, a complete graph) every vertex has the
same eccentricity, so the periphery is the whole vertex set. The periphery is
`{}` unless the graph is connected (strongly connected, if directed).
