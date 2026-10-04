# UndirectedGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`UndirectedGraph[g] gives the undirected graph underlying g: u->v and v->u become the single edge u<->v, whose weight is the sum of theirs. Edges are oriented and ordered by VertexList position. An undirected g is returned unchanged.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= EdgeList[UndirectedGraph[Graph[{1->2,2->1,2->3}]]]
Out[1]= {1 <-> 2, 2 <-> 3}

In[2]:= EdgeList[UndirectedGraph[Graph[{3->1, 2->3}]]]
Out[2]= {3 <-> 1, 3 <-> 2}
```

### Options (1)

```mathematica
In[3]:= InputForm[UndirectedGraph[Graph[{1,2,3},{1->2,2->1,3->1},EdgeWeight->{2,3,4}]]]
Out[3]= Graph[{1, 2, 3}, {1 <-> 2, 1 <-> 3}, EdgeWeight -> {5, 4}]
```

### Applications (3)

A directed triangle becomes undirected

```mathematica
In[4]:= EdgeList[UndirectedGraph[Graph[{1, 2, 3}, {1 -> 2, 2 -> 3, 3 -> 1}]]]
Out[4]= {1 <-> 2, 1 <-> 3, 2 <-> 3}
```

Direction dropped

```mathematica
In[5]:= UndirectedGraphQ[UndirectedGraph[Graph[{1, 2}, {1 -> 2}]]]
Out[5]= True
```

The two arcs merge into one edge

```mathematica
In[6]:= EdgeList[UndirectedGraph[Graph[{1, 2}, {1 -> 2, 2 -> 1}]]]
Out[6]= {1 <-> 2}
```

## Implementation notes

**Algorithm.** `builtin_undirected_graph` gives the undirected graph underlying
`g`: each edge loses its direction and the opposite arcs `u -> v` and `v -> u`
merge into a single `u <-> v` whose weight is the **sum** of theirs. An undirected
`g` is returned unchanged. The builder maps each edge to a normalized `(lo, hi)`
endpoint pair (lower `VertexList` position first), stably counting-sorts the
edges by `(lo, hi)` (row-major, upper triangle — Mathematica's order), then walks
runs of equal pairs: a run of one keeps its weight, a longer run has its weights
summed by building and evaluating a `Plus`. Because summing weights can run
arbitrary code and evict `g` from the memo, it first takes private endpoint
copies (`gops_view_own`).

**Data structures.** A `GopsView` (owned when weighted), per-edge `perm`/`lo`/`hi`
arrays and the subsystem counting sort (`gops_csort`), fresh vertex/edge/weight/
endpoint arrays, and a shared interned `UndirectedEdge` head. `gops_graph_new`
seeds the result.

**Complexity / limits.** `O(V + E)` apart from the stable sort. The opposite
form, `DirectedGraph`, replaces each `u <-> v` by the arc pair `u -> v`, `v ->
u`.

- `Protected`. A non-graph argument is left unevaluated.
- `u -> v` and `v -> u` merge into one edge whose weight is the sum (`Plus`) of theirs.
- Edges are oriented and ordered by `VertexList` position (upper triangle, row-major).
- An undirected `g` is returned unchanged.
- Performance: 1.1–1.6x faster than Mathematica 15 at `10^5` vertices
  (`benchmarks/93-graph-ops-editing`).

**Attributes:** `Protected`.

## References

**See also:** [Plus](../../arithmetic/Plus/), [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/gops_transform.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_transform.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`UndirectedGraph[g]` drops direction from every edge. A pair of opposite arcs
`u -> v` and `v -> u` collapses into the single undirected edge `u <-> v`; when
the graph is weighted, the merged edge's weight is the sum of the two arcs'
weights.

Edges come out oriented and ordered by `VertexList` position (upper triangle,
row-major). An already-undirected graph is returned unchanged. The inverse-
direction operation is `DirectedGraph`.
