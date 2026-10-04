# UndirectedGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`UndirectedGraphQ[g] gives True if every edge of g is undirected (including when g has no edges), and False otherwise.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= UndirectedGraphQ[Graph[{1,2},{}]]
Out[1]= True

In[2]:= UndirectedGraphQ[CycleGraph[3]]
Out[2]= True

In[3]:= UndirectedGraphQ[Graph[{1,2,3},{1->2,2<->3}]]
Out[3]= False

In[4]:= UndirectedGraphQ[x]
Out[4]= False
```

### Applications (5)

Every edge of a cycle graph is undirected

```mathematica
In[5]:= UndirectedGraphQ[CycleGraph[5]]
Out[5]= True
```

A single directed edge

```mathematica
In[6]:= UndirectedGraphQ[Graph[{1 -> 2}]]
Out[6]= False
```

Explicitly undirected edges

```mathematica
In[7]:= UndirectedGraphQ[Graph[{1 <-> 2, 2 <-> 3}]]
Out[7]= True
```

The generators build undirected graphs

```mathematica
In[8]:= UndirectedGraphQ[CompleteGraph[4]]
Out[8]= True
```

A non-graph argument is False

```mathematica
In[9]:= UndirectedGraphQ[7]
Out[9]= False
```

## Implementation notes

**Algorithm.** `builtin_undirected_graph_q` returns `True` exactly when `g` has no
directed edges: it is the single test `graph_directed_edge_count(g) == 0`. An edgeless
graph is therefore undirected (there is nothing oriented), matching the Wolfram Language.
A graph that mixes directed and undirected edges, or is purely directed, gives `False`.

**Data structures.** `graph_directed_edge_count` is an `O(1)` query against `g`'s
validated-graph memo (the per-node cache described in `graph.h`): on the first call it
walks `g`'s edge list once to count `DirectedEdge` nodes and memoizes the vertex index and
edge-key set; later calls read the cached count directly. The graph is the plain
`Graph[List verts, List edges]` `Expr` tree — no separate representation is materialized.

**Complexity / limits.** `O(1)` on a memo hit, one `O(E)` pass on the first query. A
non-graph argument makes `graph_directed_edge_count` return `−1`, so the `== 0` test is
`False` — the predicate never leaves the call unevaluated.

- `Protected`. Edgeless graphs are undirected; a mixed graph is neither directed
  nor undirected (see `DirectedGraphQ`).
- Shared by all the structural predicates (`UndirectedGraphQ`, `EmptyGraphQ`,
  `CompleteGraphQ`, `BipartiteGraphQ`, `VertexQ`, `EdgeQ`, `AcyclicGraphQ`,
  `TreeGraphQ`): each gives `False` — never unevaluated — for an argument that
  is not a valid graph, and each runs in linear time on first use and `O(1)`
  when repeated on the same graph (see the *Performance model* section of this
  file's preamble).

**Attributes:** `Protected`.

## References

**See also:** [DirectedGraphQ](../../graphs/DirectedGraphQ/), [EmptyGraphQ](../../graphs/EmptyGraphQ/), [CompleteGraphQ](../../graphs/CompleteGraphQ/), [BipartiteGraphQ](../../graphs/BipartiteGraphQ/), [VertexQ](../../graphs/VertexQ/), [EdgeQ](../../graphs/EdgeQ/), [AcyclicGraphQ](../../graphs/AcyclicGraphQ/), [TreeGraphQ](../../graphs/TreeGraphQ/)

- Source: [`src/graph/graphprops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graphprops.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

`UndirectedGraphQ[g]` is `True` precisely when `g` has no directed edges, which includes
the edgeless case: a graph with no edges is undirected because nothing in it is oriented.
A purely directed graph, or one that mixes the two edge kinds, gives `False`.

The test is a single memoized directed-edge count, so it is `O(1)` after the first query
on a given graph. A non-graph argument yields `False` rather than remaining unevaluated.
