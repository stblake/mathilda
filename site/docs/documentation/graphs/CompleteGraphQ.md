# CompleteGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CompleteGraphQ[g] gives True if every pair of distinct vertices of g is joined by an edge in both directions (an undirected edge, or directed edges both ways). CompleteGraphQ[g, vlist] tests the subgraph induced by vlist; False if some element of vlist is not a vertex of g.`**

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= CompleteGraphQ[Graph[{1->2,2->1}]]
Out[1]= True

In[2]:= CompleteGraphQ[Graph[{1->2}]]
Out[2]= False

In[3]:= CompleteGraphQ[CompleteGraph[5]]
Out[3]= True

In[4]:= CompleteGraphQ[CycleGraph[4], {1,2}]
Out[4]= True

In[5]:= CompleteGraphQ[CycleGraph[4], {1,3}]
Out[5]= False

In[6]:= CompleteGraphQ[CycleGraph[4], {1,9}]
Out[6]= False
```

### Applications (5)

The complete graph on 4 vertices

```mathematica
In[7]:= CompleteGraphQ[CompleteGraph[4]]
Out[7]= True
```

A 5-cycle is far from complete

```mathematica
In[8]:= CompleteGraphQ[CycleGraph[5]]
Out[8]= False
```

The subgraph induced by three vertices is still complete

```mathematica
In[9]:= CompleteGraphQ[CompleteGraph[4], {1, 2, 3}]
Out[9]= True
```

Those three vertices are not pairwise adjacent

```mathematica
In[10]:= CompleteGraphQ[CycleGraph[5], {1, 2, 3}]
Out[10]= False
```

A non-graph argument is simply False

```mathematica
In[11]:= CompleteGraphQ[5]
Out[11]= False
```

## Implementation notes

**Algorithm.** `builtin_complete_graph_q` tests whether every ordered pair of
distinct vertices `(u, v)` is joined by an edge usable from `u` to `v` — an
`UndirectedEdge` between them, or a `DirectedEdge[u, v]`, so a complete *directed*
graph needs both `u -> v` and `v -> u`. It first reads `graph_directed_edge_count`
(an `O(1)` memo query that also validates `g`); a non-graph gives that count `< 0`
and the head returns `False`. For a simple graph with a single edge kind, completeness
reduces to an edge count — `n(n−1)/2` undirected edges or `n(n−1)` directed ones, since
parallel edges cannot occur — so no traversal is needed. A mixed graph (both edge kinds)
may pair `u<->v` with `u->v`, so it falls back to the adjacency scan `induced_complete`,
which stamps the distinct out-neighbours of each vertex and checks that each reaches the
other `k−1`. Graphs on 0 or 1 vertices are complete (no pair can fail).

**Data structures.** `CompleteGraphQ[g, vlist]` builds the integer-indexed `GraphAdj`
(successor/predecessor CSR) and a `char` selection mask `sel[]` over the vertices named
in `vlist`, resolving each through `graph_vertex_position` (the memoized vertex index);
if any element of `vlist` is not a vertex of `g`, the result is `False`. `induced_complete`
uses a per-`u` `stamp[]` array so an undirected edge listed alongside a redundant directed
one is not double-counted. The graph itself is the ordinary `Graph[List verts, List edges]`
`Expr` tree; validation and the edge count come from its per-node memo.

**Complexity / limits.** The single-edge-kind fast path is `O(1)` after validation; the
induced/mixed scan is `O(V + E)`. An allocation failure in the scan leaves the call
unevaluated (`NULL`) rather than answering `False`. Every non-graph argument returns
`False`, never an unevaluated expression.

- `Protected`. A pair `(u, v)` is covered by an undirected edge or a directed
  `u -> v`, so a complete directed graph needs both directions.
- Graphs with 0 or 1 vertices are complete.
- `CompleteGraphQ[g, vlist]` gives `False` if some element is not a vertex of
  `g`; repeats are ignored; `{}` is complete.
- `False` for a non-graph (see `UndirectedGraphQ`).

**Attributes:** `Protected`.

## References

**See also:** [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- Source: [`src/graph/graphprops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graphprops.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

Completeness means every pair of distinct vertices is joined; for a directed graph that
requires edges in *both* directions (equivalently, an undirected edge). Graphs on zero or
one vertices are complete by default, since there is no pair that could fail.

The two-argument form `CompleteGraphQ[g, vlist]` tests the subgraph of `g` induced by
`vlist`, and gives `False` if any element of `vlist` is not actually a vertex of `g`.
Like every `*Q` predicate, a non-graph argument yields `False` rather than staying
unevaluated.
