# GraphRadius

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphRadius[g] gives the minimum vertex eccentricity of g. Infinity unless g is (strongly) connected.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= GraphDiameter[PathGraph[5]]
Out[1]= 4

In[2]:= GraphRadius[PathGraph[5]]
Out[2]= 2

In[3]:= GraphDiameter[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= Infinity

In[4]:= GraphDiameter[Graph[{},{}]]
Out[4]= 0
```

### Options (2)

```mathematica
In[5]:= GraphRadius[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}]]
Out[5]= 3.0

In[6]:= GraphDiameter[Graph[{1,2,3},{1<->2,2<->3,1<->3}, EdgeWeight->{1,1,5}]]
Out[6]= 2.0
```

### Applications (3)

The central vertex has eccentricity 2

```mathematica
In[7]:= GraphRadius[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 5}]]
Out[7]= 2
```

```mathematica
In[8]:= gw = Graph[{1 <-> 2, 2 <-> 3}, EdgeWeight -> {1, 4}];
```

Weighted: a machine real

```mathematica
In[9]:= GraphRadius[gw]
Out[9]= 4.0
```

## Implementation notes

**Algorithm.** `builtin_graph_radius` is the minimum vertex eccentricity. Through the shared
`extremal_compute`, it reads a cached per-source distance summary (`gmet_dist_summary`:
multi-source BFS when unweighted, Dijkstra when weighted), computes each vertex's eccentricity
`ecc[i]` (the largest distance from `i` to a vertex it reaches), and returns the minimum. An
`O(V+E)` strong-connectivity test short-circuits first: an unweighted graph that is not
(strongly) connected has radius `Infinity`. A weighted graph can still return a finite radius
even when some pair is unreachable, following Wolfram's rule that a weighted eccentricity is
`Infinity` only when the vertex fails to reach all `n-1` others.

**Data structures.** The memoised distance summary's `ecc[]`/`reach[]` arrays, shared with
`GraphCenter` and `MeanGraphDistance` so the three heads cost one traversal between them. The
value is returned exactly (an `Integer`) for an unweighted graph and as a machine `Real` for a
weighted one, through `gmet_distance_value`. A weighted eccentricity tie uses a relative
tolerance `ECC_TIE_TOL = 1e-12`.

**Complexity / limits.** `O(V(V+E))` unweighted, `O(V(E + V log V))` weighted. The empty graph
returns `0`; a non-(strongly-)connected unweighted graph returns `Infinity`; a non-graph
argument returns unevaluated.

- *(w)* weight-aware (machine reals when weighted).
- Unweighted and not strongly connected (not connected, if undirected):
  `Infinity`.
- Weighted: taken over the weighted eccentricities, so a weighted digraph that
  is not strongly connected can have a finite radius
  (*reverse-engineered*).
- Graphs with no vertices give diameter/radius `0`.
- Reduced from the cached MS-BFS per-source summary, so calling several of
  `GraphDiameter`, `GraphRadius`, `GraphCenter`, `GraphPeriphery`,
  `ClosenessCentrality`, `EccentricityCentrality` on one graph pays for one
  all-pairs pass.

**Attributes:** `Protected`.

## References

**See also:** [GraphDiameter](../../graphs/GraphDiameter/), [GraphCenter](../../graphs/GraphCenter/), [GraphPeriphery](../../graphs/GraphPeriphery/), [ClosenessCentrality](../../graphs/ClosenessCentrality/), [EccentricityCentrality](../../graphs/EccentricityCentrality/)

- Source: [`src/graph/gmet_distance.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_distance.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The radius is the minimum eccentricity over all vertices (the eccentricity of a vertex being its
greatest distance to any other). It is an exact `Integer` for an unweighted graph and a machine
`Real` for a weighted one.

It is `Infinity` unless the graph is connected (strongly connected, for a directed graph).
`GraphCenter` returns the vertices that attain this minimum eccentricity.
