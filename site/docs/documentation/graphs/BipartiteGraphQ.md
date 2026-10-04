# BipartiteGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BipartiteGraphQ[g] gives True if the vertices of g split into two sets with every edge running between them (edge direction is ignored).`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= BipartiteGraphQ[CycleGraph[5]]
Out[1]= False

In[2]:= BipartiteGraphQ[CycleGraph[6]]
Out[2]= True

In[3]:= BipartiteGraphQ[Graph[{1,2,3},{}]]
Out[3]= True

In[4]:= BipartiteGraphQ[x]
Out[4]= False
```

### Applications (4)

An even cycle is bipartite

```mathematica
In[5]:= BipartiteGraphQ[CycleGraph[4]]
Out[5]= True
```

An odd cycle is not

```mathematica
In[6]:= BipartiteGraphQ[CycleGraph[5]]
Out[6]= False
```

Every complete bipartite graph qualifies

```mathematica
In[7]:= BipartiteGraphQ[CompleteGraph[{2, 3}]]
Out[7]= True
```

K4 contains a triangle, so False

```mathematica
In[8]:= BipartiteGraphQ[CompleteGraph[4]]
Out[8]= False
```

## Implementation notes

**Algorithm.** `builtin_bipartite_graph_q` tests whether the *underlying
undirected* graph is 2-colourable. It runs a breadth-first 2-colouring over
every component: each vertex is given a side (0/1) and a neighbour on the same
side is a conflict, which proves an odd cycle and so non-bipartiteness. Edge
direction is ignored — the neighbourhood of a vertex is both its successors
(`out[]`) and its predecessors (`in[]`) — because bipartiteness is a property of
the vertex split, not of orientation. An edgeless graph is trivially bipartite;
a non-graph argument gives `False` (never unevaluated), as every `*Q` predicate
does.

**Data structures.** It builds the integer-indexed `GraphAdj` (CSR successor and
predecessor lists, `src/graph/graph.h`), a `signed char side[]` array, and an
explicit BFS queue `q[]`. The result is memoized on the graph node through
`graph_prop_set(g, GRAPH_PROP_BIPARTITE, ...)`, so a repeat query on the same
`Graph` object is `O(1)` — the same property-caching an atomic Mathematica
`Graph` does.

**Complexity / limits.** `O(V + E)` for the first call, `O(1)` thereafter.
Allocation failure is the only path that returns `NULL`; a malformed graph
returns `False`.

- `Protected`. Edge direction is ignored; edgeless graphs are bipartite.
  `False` for a non-graph (see `UndirectedGraphQ`).

**Attributes:** `Protected`.

## References

**See also:** [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- Source: [`src/graph/graphprops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graphprops.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

A graph is bipartite exactly when it has no odd cycle. The test is a BFS
2-colouring over each connected component; direction is ignored, so a directed
graph is tested as its underlying undirected graph.

The answer is cached on the graph object, so repeated bipartiteness queries on
the same `Graph` cost nothing after the first. A non-graph argument gives
`False` rather than staying unevaluated.
