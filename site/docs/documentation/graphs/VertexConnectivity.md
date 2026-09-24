# VertexConnectivity

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexConnectivity[g] gives the minimum number of vertices whose removal disconnects g.`**

## Examples (5)

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

## Algorithm

connectivity.c - ConnectedGraphQ[g] and VertexConnectivity[g].

Both operate on the underlying undirected graph.

```text
  ConnectedGraphQ[g]    True iff g has >= 1 vertex and forms a single
                        connected component.
  VertexConnectivity[g] the least number of vertices whose removal
                        disconnects g (n-1 for a complete graph, 0 if already
                        disconnected or trivial). Computed by brute-force
                        search over vertex subsets -- exact, intended for the
                        small graphs of a pico-CAS.
```

Memory (SPEC section 4): returns freshly-allocated results; frees res.

## Implementation notes

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

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
