# VertexConnectivity

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexConnectivity[g] gives the minimum number of vertices whose removal disconnects g.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= VertexConnectivity[CycleGraph[5]]
Out[1]= 2

In[2]:= VertexConnectivity[CompleteGraph[5]]
Out[2]= 4

In[3]:= VertexConnectivity[StarGraph[5]]
Out[3]= 1

In[4]:= VertexConnectivity[Graph[{1,2,3,4},{1<->2,3<->4}]]
Out[4]= 0

In[5]:= VertexConnectivity[Graph[{1->2,2->3,3->4,4->1}]]
Out[5]= 2
```

### Applications (3)

A complete graph on n needs n-1 removed

```mathematica
In[6]:= VertexConnectivity[CompleteGraph[4]]
Out[6]= 3
```

Removing one interior vertex breaks a path

```mathematica
In[7]:= VertexConnectivity[PathGraph[5]]
Out[7]= 1
```

Vertex 4 is a cut vertex

```mathematica
In[8]:= VertexConnectivity[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 2 <-> 4, 4 <-> 5}]]
Out[8]= 1
```

## Algorithm

connectivity.c - ConnectedGraphQ[g] and VertexConnectivity[g].

```text
  ConnectedGraphQ[g]    True iff g has >= 1 vertex and forms a single
                        connected component: strongly connected when g has
                        a directed edge (Mathematica's rule), connected
                        otherwise.
  VertexConnectivity[g] the least number of vertices whose removal
                        disconnects g (n-1 for a complete graph, 0 if already
                        disconnected or trivial). Computed by brute-force
                        search over vertex subsets -- exact, intended for the
                        small graphs of a pico-CAS.
```

Memory (SPEC section 4): returns freshly-allocated results; frees res.

## Implementation notes

**Algorithm.** `builtin_vertex_connectivity` gives the minimum number of vertices whose removal
disconnects the graph. The primary engine (`galg_vertex_connectivity` in `galg_flow.c`) is a
max-flow method: Even's split-vertex construction turns each vertex `v` into `v_in -> v_out`
with capacity 1, so a vertex cut becomes a minimum `s`-`t` cut, and Esfahanian–Hakimi pair
selection bounds the number of local `kappa(s,t)` computations needed for the global minimum.
Each local value is a Dinic max flow on the unit-capacity residual network. A brute-force
`C(n,k)` subset enumeration (removing every `k`-subset and testing connectivity) is retained
only as an allocation-failure fallback.

**Data structures.** The split residual network as a CSR of unit-capacity arcs, with Dinic's
BFS level array and current-arc pointers. The fallback uses a `GraphAdj`, a combination index
array, and a `removed[]` mask passed to `graph_count_components`.

**Complexity / limits.** The flow path is polynomial (each Dinic max flow on unit capacities is
about `O(E sqrt(V))`); the fallback is exponential. A complete graph returns `n-1`; an
already-disconnected graph or one with `n <= 1` returns `0`. A non-graph argument returns
unevaluated.

- `Protected`. Gives `n-1` for `K_n`, and `0` if `g` is already disconnected.
- Algorithm: originally an exact brute-force search over vertex subsets,
  intended for small graphs. It now uses the max-flow machinery shared with
  `FindVertexCut` (via `galg_vertex_connectivity`: Even's split-vertex network
  with Esfahanian-Hakimi pair selection), replacing the former exponential
  subset search.
- Edge direction is currently ignored (the underlying undirected graph is
  used), so a directed cycle scores like an undirected one. Mathematica uses
  strong connectivity for directed graphs (a directed cycle has vertex
  connectivity 1).

**Attributes:** `Protected`.

## References

**See also:** [FindVertexCut](../../graphs/FindVertexCut/)

- S. Even, *An algorithm for determining whether the connectivity of a graph is at least k*, SIAM J. Comput. **4** (1975) 393-396.
- A. H. Esfahanian and S. L. Hakimi, *On computing the connectivities of graphs and digraphs*, Networks **14** (1984) 355-366.
- Source: [`src/graph/connectivity.c`](https://github.com/stblake/mathilda/blob/main/src/graph/connectivity.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

The result is the smallest number of vertices whose removal disconnects the graph (or reduces it
to a single vertex). A complete graph on `n` vertices has connectivity `n-1`; a graph that is
already disconnected, or has at most one vertex, has connectivity `0`.

`FindVertexCut` returns an actual minimum separating set of that size.
