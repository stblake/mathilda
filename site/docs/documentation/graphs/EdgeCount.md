# EdgeCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeCount[g] gives the number of edges in the graph g.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= VertexCount[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]
Out[1]= 4

In[2]:= EdgeCount[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]
Out[2]= 4

In[3]:= VertexCount[Graph[{}, {}]]
Out[3]= 0

In[4]:= EdgeCount[Graph[{1,2,3},{}]]
Out[4]= 0

In[5]:= EdgeCount[x]
Out[5]= EdgeCount[x]
```

### Applications (3)

A complete graph on 4 vertices has binomial(4,2) edges

```mathematica
In[6]:= EdgeCount[CompleteGraph[4]]
Out[6]= 6
```

A cycle has as many edges as vertices

```mathematica
In[7]:= EdgeCount[CycleGraph[5]]
Out[7]= 5
```

Directed edges count once each

```mathematica
In[8]:= EdgeCount[Graph[{1 -> 2, 2 -> 3}]]
Out[8]= 2
```

## Implementation notes

**Algorithm.** `builtin_edge_count` is a thin reader over the canonical form
`Graph[List[verts], List[edges]]`: it validates `g` and returns the integer
`arg_count` of the edge sublist `args[1]`. Every edge in the canonical list —
whether a `DirectedEdge` or an `UndirectedEdge` — is counted exactly once, so
the count is the length of `EdgeList[g]`. When the argument is not a valid graph
the call is forwarded to `hyp_edge_count`, which answers for a `Hypergraph` and
otherwise returns `NULL` (the call stays unevaluated).

**Data structures.** None beyond the result. The edge list is already a stored
child of the graph node, so no traversal or allocation happens past constructing
the fresh integer the evaluator then owns.

**Complexity / limits.** `O(1)` — the arity is read straight from the stored
list header. The one-time validation `graph_is_valid(g)` is `O(1)` on a memo hit
and `O(V + E)` the first time a graph is seen. This head takes exactly one
argument; any other arity returns `NULL`, so there is no edge-pattern counting
form (`EdgeCount[g, patt]` is left unevaluated).

- `Protected`. Cardinalities read from the canonical form; unevaluated on a
  non-graph (see `VertexList`).

**Attributes:** `Protected`.

## References

**See also:** [VertexCount](../../graphs/VertexCount/), [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/counts.c`](https://github.com/stblake/mathilda/blob/main/src/graph/counts.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`EdgeCount[g]` is the length of `EdgeList[g]`: every edge of the canonical form,
directed or undirected, is counted exactly once.

The head takes a single argument and declines otherwise — there is no
edge-pattern counting form, so `EdgeCount[g, patt]` is left unevaluated. A
non-graph argument that is a `Hypergraph` is answered with its hyperedge count;
anything else stays unevaluated.
