# FindVertexCut

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindVertexCut[g] gives a minimum set of vertices whose removal disconnects the underlying undirected graph of g (n-1 vertices for a complete graph); FindVertexCut[g, s, t] gives a minimum set separating s from t ({} when they are adjacent).`**

## Examples (6)

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

## Implementation notes

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

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
