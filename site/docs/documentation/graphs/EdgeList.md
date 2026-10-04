# EdgeList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeList[g] gives the list of edges of the graph g.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EdgeList[Graph[{1,2,3},{1->2, 2<->3}]]
Out[1]= {1 -> 2, 2 <-> 3}

In[2]:= EdgeList[CycleGraph[4]]
Out[2]= {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}

In[3]:= EdgeList[Graph[{1,2},{}]]
Out[3]= {}

In[4]:= EdgeList[5]
Out[4]= EdgeList[5]
```

### Applications (3)

Undirected edges, in canonical order

```mathematica
In[5]:= EdgeList[CompleteGraph[3]]
Out[5]= {1 <-> 2, 1 <-> 3, 2 <-> 3}
```

Directed edges keep their arrows

```mathematica
In[6]:= EdgeList[Graph[{1 -> 2, 2 -> 3}]]
Out[6]= {1 -> 2, 2 -> 3}
```

The count is the list length

```mathematica
In[7]:= EdgeCount[CycleGraph[4]] == Length[EdgeList[CycleGraph[4]]]
Out[7]= True
```

## Implementation notes

**Algorithm.** `builtin_edge_list` validates `g` and returns a deep copy of the
graph's canonical edge sublist `args[1]`, preserving stored order. Each element
is already in the canonical two-argument `DirectedEdge[u, v]` / `UndirectedEdge[u, v]`
form the constructor normalises edge sugar (`u -> v`, `u <-> v`) into, so
`EdgeList` is the inverse of the edge half of `Graph[...]`. When the argument is
not a valid graph the call is forwarded to `hyp_edge_list` (which returns the
hyperedges of a `Hypergraph`) and otherwise yields `NULL`.

**Data structures.** The result is `expr_copy(g->data.function.args[1])` — an
independent `List` of edge nodes, so the caller can mutate or consume it without
disturbing the graph. No adjacency structure is consulted.

**Complexity / limits.** `O(E)` for the copy, after an `O(1)` memo-hit (or
`O(V + E)` first-time) validation. One-argument only; other arities return
`NULL`.

- `Protected`. A thin reader over the canonical form; unevaluated on a non-graph
  (see `VertexList`).
- Edges print in operator form (`1 -> 2`, `2 <-> 3`), but are
  `DirectedEdge`/`UndirectedEdge` internally.

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/edgelist.c`](https://github.com/stblake/mathilda/blob/main/src/graph/edgelist.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`EdgeList[g]` returns the graph's edges in the canonical two-argument form:
`UndirectedEdge[u, v]` (printed `u <-> v`) and `DirectedEdge[u, v]` (printed
`u -> v`). Edge sugar given to the constructor is already normalised to these
heads, so `EdgeList` round-trips the edge half of `Graph[...]`.

A `Hypergraph` argument returns its hyperedges; any other non-graph argument
leaves the call unevaluated.
