# GraphDiameter

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphDiameter[g] gives the greatest distance between two vertices of g: the maximum vertex eccentricity. Infinity unless g is (strongly) connected.`**

## Examples (11)

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

### Applications (5)

The farthest pair on a six-cycle is three steps apart

```mathematica
In[7]:= GraphDiameter[CycleGraph[6]]
Out[7]= 3
```

A path of five vertices spans four edges

```mathematica
In[8]:= GraphDiameter[PathGraph[{1, 2, 3, 4, 5}]]
Out[8]= 4
```

Corner to corner is (3-1) + (4-1) steps

```mathematica
In[9]:= GraphDiameter[GridGraph[{3, 4}]]
Out[9]= 5
```

An isolated vertex makes the graph disconnected

```mathematica
In[10]:= GraphDiameter[Graph[{1, 2, 3}, {UndirectedEdge[1, 2]}]]
Out[10]= Infinity
```

Distances follow edge direction on a directed cycle

```mathematica
In[11]:= GraphDiameter[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[2, 3], DirectedEdge[3, 1]}]]
Out[11]= 2
```

## Implementation notes

**Algorithm.** `builtin_graph_diameter` is the `EX_DIAMETER` case of the shared `extremal` routine. It validates the graph and reads any `EdgeWeight` list, then runs an `O(n + m)` strong-connectivity test over the out-arcs. An unweighted graph that is not strongly connected (connected, when undirected) answers `Infinity` immediately, with no all-pairs work. Otherwise it reads the cached per-source summary (`gmet_dist_summary`) and returns the largest eccentricity. Unweighted distances come from a bit-parallel multi-source BFS, 256 sources per adjacency sweep; weighted graphs use a binary-heap Dijkstra per source and a machine-real result. The empty graph gives `0`.

**Data structures.** The graph is the Expr tree `Graph[List[v1, ...], List[edge1, ...]]`; `gmet_csr_build` turns its memoized edge-index view (`graph_edge_indices`) into a CSR arc list, with an undirected edge stored as two arcs. The summary holds per-source eccentricity, reach count and distance sum arrays, so `GraphRadius`, `GraphCenter` and `MeanGraphDistance` reuse the same traversal. Results are cached by expression.

**Complexity / limits.** Unweighted cost is `O(n (n + m) / 256)` word operations, with source batches spread over the thread team; weighted cost is `O(n (m + n) log n)`. A disconnected unweighted graph is answered in `O(n + m)`. Weighted graphs return `Infinity` whenever some pair is unreachable. Symbolic, complex or negative weights leave the call unevaluated.

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

**See also:** [GraphRadius](../../graphs/GraphRadius/), [GraphCenter](../../graphs/GraphCenter/), [GraphPeriphery](../../graphs/GraphPeriphery/), [ClosenessCentrality](../../graphs/ClosenessCentrality/), [EccentricityCentrality](../../graphs/EccentricityCentrality/)

- Source: [`src/graph/gmet_distance.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_distance.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The diameter is the largest shortest-path distance over all ordered pairs of vertices. An unweighted graph gives an exact integer; a graph carrying `EdgeWeight` gives a machine real.

If some vertex cannot be reached from another, the diameter is `Infinity`. For a directed graph this means the graph is not strongly connected, which is checked first in linear time so the disconnected case costs no all-pairs work.
