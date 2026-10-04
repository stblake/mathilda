# FindVertexCut

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindVertexCut[g] gives a minimum set of vertices whose removal disconnects the underlying undirected graph of g (n-1 vertices for a complete graph); FindVertexCut[g, s, t] gives a minimum set separating s from t ({} when they are adjacent).`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= FindVertexCut[CycleGraph[6]]
Out[1]= {2, 4}

In[2]:= FindVertexCut[CycleGraph[6], 1, 4]
Out[2]= {3, 5}

In[3]:= FindVertexCut[CycleGraph[6], 1, 2]
Out[3]= {}

In[4]:= FindVertexCut[CompleteGraph[4]]
Out[4]= {1, 2, 3}

In[5]:= FindVertexCut[Graph[{1->2,UndirectedEdge[1,3],UndirectedEdge[2,3]}]]
Out[5]= {}

In[6]:= VertexConnectivity[CycleGraph[6]]
Out[6]= 2
```

### Applications (2)

One cut vertex suffices

```mathematica
In[7]:= FindVertexCut[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 2 <-> 4, 4 <-> 5}]]
Out[7]= {2}
```

A complete graph needs n-1 vertices removed

```mathematica
In[8]:= FindVertexCut[CompleteGraph[4]]
Out[8]= {1, 2, 3}
```

## Implementation notes

**Algorithm.** `builtin_find_vertex_cut` returns an actual minimum set of vertices whose removal
disconnects the underlying undirected graph (or, in the three-argument form, separates `s` from
`t`). It uses the same Even split-vertex construction and Esfahanian–Hakimi pair selection as
`VertexConnectivity`, but then extracts the separator from the Dinic max flow: the cut vertices
are those whose in-copy can no longer reach `t_in` while their out-copy still can. Ties are
broken toward the cut closest to `t`.

**Data structures.** A sorted-CSR undirected view `GalgUG` of the graph, the split residual
network `GfNet`, a `sep[]` output index array, and a `mark[]` reachability array over the
residual graph after the final flow.

**Complexity / limits.** Polynomial — Dinic max flow on unit capacities, as for
`VertexConnectivity`. `FindVertexCut[g]` on a complete undirected graph returns its first `n-1`
vertices; `FindVertexCut[g, s, t]` with `s` and `t` adjacent returns `{}` (they cannot be
separated without deleting an endpoint). A non-graph argument, or `s == t`, returns unevaluated.

- `Protected`; unevaluated on a non-graph.
- Works on the underlying undirected graph; vertices are returned in
  `VertexList` order.
- The `s`-`t` separator is the one closest to `t`; `{}` for adjacent `s`, `t`.
- Mathematica's conventions: a complete undirected graph gives its first `n-1`
  vertices; a graph with a directed edge whose underlying graph is complete
  gives `{}`.
- Algorithm: Even's split-vertex network with Esfahanian-Hakimi pair selection,
  solved by Dinic max flow. `VertexConnectivity` uses the same machinery (via
  `galg_vertex_connectivity`), replacing its former exponential subset search
  (identical answers on 300 random graphs).

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/), [VertexConnectivity](../../graphs/VertexConnectivity/)

- S. Even, *An algorithm for determining whether the connectivity of a graph is at least k*, SIAM J. Comput. **4** (1975) 393-396.
- E. A. Dinic, *Algorithm for solution of a problem of maximum flow in networks with power estimation*, Soviet Math. Dokl. **11** (1970) 1277-1280.
- Source: [`src/graph/galg_flow.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_flow.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

`FindVertexCut[g]` gives a smallest vertex set whose removal disconnects the underlying
undirected graph; its length equals `VertexConnectivity[g]`. For a complete graph on `n`
vertices this is `n-1` vertices (you cannot disconnect a clique without deleting all but one).

`FindVertexCut[g, s, t]` gives a smallest set separating `s` from `t`, and returns `{}` when `s`
and `t` are adjacent.
