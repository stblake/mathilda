# EccentricityCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EccentricityCentrality[g] gives 1/VertexEccentricity[g, v] for each vertex v (0 when the eccentricity is 0). Uses EdgeWeight as lengths.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= EccentricityCentrality[PathGraph[5]]
Out[1]= {0.25, 0.333333, 0.5, 0.333333, 0.25}

In[2]:= EccentricityCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[2]= {0.333333, 0.5, 0.5, 0.0}

In[3]:= EccentricityCentrality[Graph[{1,2},{}]]
Out[3]= {0.0, 0.0}
```

### Options (1)

```mathematica
In[4]:= EccentricityCentrality[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}]]
Out[4]= {0.333333, 0.5, 0.0}
```

### Applications (4)

Every vertex of a cycle has the same eccentricity

```mathematica
In[5]:= EccentricityCentrality[CycleGraph[5]]
Out[5]= {0.5, 0.5, 0.5, 0.5, 0.5}
```

The middle vertex is the most central

```mathematica
In[6]:= EccentricityCentrality[PathGraph[{1, 2, 3, 4, 5}]]
Out[6]= {0.25, 0.333333, 0.5, 0.333333, 0.25}
```

The hub reaches every leaf in one step

```mathematica
In[7]:= EccentricityCentrality[StarGraph[5]]
Out[7]= {1.0, 0.5, 0.5, 0.5, 0.5}
```

Edge weights are used as path lengths

```mathematica
In[8]:= EccentricityCentrality[Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight -> {1, 2}]]
Out[8]= {0.333333, 0.5, 0.333333}
```

## Implementation notes

**Algorithm.** `builtin_eccentricity_centrality` is a one-line call to
`summary_centrality(res, 1)`. That helper fetches the per-source distance summary
`gmet_dist_summary(g)` and returns, for every vertex `v`, the reciprocal
eccentricity `1/ecc[v]` (and `0` when `ecc[v] = 0`, i.e. an isolated vertex), where
`ecc[v]` is the largest shortest-path distance from `v` to any vertex reachable
along its out-arcs. The summary itself does all the work: it builds an out-arc CSR
and, for each source, records the number of reachable vertices, the sum of their
distances, and their maximum. Unweighted graphs use a bit-parallel multi-source
BFS (a block of sources advanced together with a bitset frontier, parallelised over
blocks); when an `EdgeWeight` is present it switches to a binary-heap Dijkstra per
source, using the weights as lengths.

**Data structures.** A `GmetCSR` compressed-sparse-row adjacency, and a
`GmetDistSummary` (`reach[]`, `sum[]`, `ecc[]`) cached per graph node in a small
slot table and evicted round-robin, so a second centrality query on the same graph
reuses the distances. The result is handed back through `gmet_real_vector`, which
packs into an `f64` NDArray buffer above the packing threshold and otherwise builds
a `List` of reals.

**Complexity / limits.** `O(V(V+E))` unweighted (the BFS amortised over the
word-wide source blocks) and `O(V · E log V)` weighted, both multithreaded over the
source set. The value is a machine real, never exact. A non-graph argument, or one
carrying unusable weights, leaves the call unevaluated (`NULL`). In a disconnected
graph the eccentricity is taken over the reachable set only, so a vertex's score
reflects its own component.

- *(w)* weight-aware; machine reals.
- `e(v)` is measured over the vertices reachable from `v` (also when
  weighted); the centrality is 0 when `e(v) = 0`.
- Reduced from the cached MS-BFS per-source summary (see `GraphDistanceMatrix`).

**Attributes:** `Protected`.

## References

**See also:** [GraphDistanceMatrix](../../graphs/GraphDistanceMatrix/)

- P. Hage and F. Harary, *Eccentricity and centrality in networks*, Social Networks **17** (1995) 57-63.
- Source: [`src/graph/gmet_centrality.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_centrality.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The score of a vertex `v` is `1/VertexEccentricity[g, v]`, the reciprocal of the
distance to the farthest vertex `v` can reach. A vertex at the centre of the graph
has the smallest eccentricity and therefore the largest centrality; an isolated
vertex, whose eccentricity is `0`, is assigned the score `0`.

Distances come from a BFS on an unweighted graph and from a Dijkstra shortest path
when an `EdgeWeight` is supplied, so the centrality respects edge lengths. In a
disconnected graph the eccentricity is measured over the vertices actually
reachable from `v`, i.e. within its own component.
