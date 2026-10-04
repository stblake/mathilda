# TuranGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TuranGraph[n, k] gives the Turan graph: the complete k-partite graph on n vertices with parts as equal as possible.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= TuranGraph[5, 2]
Out[1]= Graph[<5 vertices, 6 edges>]

In[2]:= EdgeList[TuranGraph[5, 2]]
Out[2]= {1 <-> 4, 1 <-> 5, 2 <-> 4, 2 <-> 5, 3 <-> 4, 3 <-> 5}

In[3]:= TuranGraph[6, 3]
Out[3]= Graph[<6 vertices, 12 edges>]

In[4]:= EdgeList[TuranGraph[4, 3]]
Out[4]= {1 <-> 3, 1 <-> 4, 2 <-> 3, 2 <-> 4, 3 <-> 4}
```

### Applications (3)

```mathematica
In[5]:= t = TuranGraph[7, 3];
```

Parts of sizes 3, 2, 2

```mathematica
In[6]:= {VertexCount[t], EdgeCount[t]}
Out[6]= {7, 16}
```

Three parts of 2 is complete tripartite K(2,2,2)

```mathematica
In[7]:= EdgeCount[TuranGraph[6, 3]]
Out[7]= 12
```

## Implementation notes

**Algorithm.** `builtin_turan_graph` builds the Turán graph `T(n, k)`: the complete
`k`-partite graph on `n` vertices with the parts as equal as possible. The part sizes are
`sz[i] = n/k + (i < n mod k ? 1 : 0)`, so the first `n mod k` parts get one extra vertex
(larger parts first). It then emits every edge between vertices in different parts via
`multipartite`. The vertices are the integers `1..n`; `k` is clamped to `n` when larger.
This is the graph that, by Turán's theorem, has the most edges among `n`-vertex graphs with no
`(k+1)`-clique.

**Data structures.** The shared generator emits edges as packed 64-bit endpoint keys into a
`Pairs` accumulator, sorts and deduplicates them, and wraps the result as an undirected
`Graph[Range[n], {...}]`. All edges share one `UndirectedEdge` head node and the integer vertex
nodes are shared by reference.

**Complexity / limits.** `O(V + E)` to emit plus `O(E log E)` to sort when not already ordered;
`E` is the Turán edge count. Guarded by a vertex cap `GEN_MAX_VERTICES = 10^8` and an edge cap
`GEN_MAX_EDGES = 5·10^7` (either returns unevaluated). Both arguments must be positive integers;
an invalid argument returns unevaluated. Always undirected.

- Undirected on `1..n`, larger parts first; the edge list is the sorted list
  of pairs `{i, j}`, `i < j`, identical to Wolfram's `EdgeList`.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- P. Turán, *On an extremal problem in graph theory* (Hungarian), Mat. Fiz. Lapok **48** (1941) 436-452.
- Source: [`src/graph/gmet_generators.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_generators.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

`TuranGraph[n, k]` is the complete `k`-partite graph on `n` vertices with parts as equal as
possible (the first `n mod k` parts get one extra vertex). By Turán's theorem it maximises the
edge count among `n`-vertex graphs that contain no `(k+1)`-clique.

The vertices are the integers `1..n` and the graph is always undirected. `TuranGraph[n, 2]` is
the balanced complete bipartite graph, and `TuranGraph[n, n]` is the complete graph.
