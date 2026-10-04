# EdgeIndex

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeIndex[g, e] gives the position of the edge e in EdgeList[g] (an undirected edge matches either orientation; u->v and u<->v sugar are accepted); EdgeIndex[g, {e1, ...}] gives a list of positions.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= VertexIndex[Graph[{a<->b, b<->c}], b]
Out[1]= 2

In[2]:= VertexIndex[Graph[{a<->b, b<->c}], {c, a}]
Out[2]= {3, 1}

In[3]:= EdgeIndex[CycleGraph[4], 3<->2]
Out[3]= 2

In[4]:= EdgeIndex[CycleGraph[4], 3->4]
Out[4]= 3

In[5]:= EdgeIndex[CycleGraph[4], {1<->2, 4<->1}]
Out[5]= {1, 4}

In[6]:= EdgeIndex[Graph[{1->2}], 2->1]
Out[6]= EdgeIndex[Graph[<2 vertices, 1 edge>], 2 -> 1]
```

### Applications (3)

Position of this edge in EdgeList

```mathematica
In[7]:= EdgeIndex[CycleGraph[4], 2 <-> 3]
Out[7]= 2
```

K3 lists {1<->2, 1<->3, 2<->3}

```mathematica
In[8]:= EdgeIndex[CompleteGraph[3], 1 <-> 3]
Out[8]= 2
```

A list of edges gives a list of positions

```mathematica
In[9]:= EdgeIndex[CycleGraph[5], {1 <-> 2, 3 <-> 4}]
Out[9]= {1, 3}
```

## Implementation notes

**Algorithm.** `builtin_edge_index` gives the 1-based position of an edge in
`EdgeList[g]`. A single edge argument — written `u -> v`, `u <-> v`, or an
explicit `DirectedEdge`/`UndirectedEdge` — is resolved to its `EdgeList`
position; an undirected edge matches either orientation. `EdgeIndex[g, {e1,
...}]` maps the resolver over a list and returns the list of positions. Any edge
not present (or a mixed-up orientation for a directed edge) leaves the call
unevaluated.

**Data structures.** It opens a `GopsView` over the graph (the shared borrowed
view of vertices, integer endpoint arrays `eu[]`/`ev[]`, and direction flags),
parses each query edge to endpoint positions with `gops_parse_edge`, and resolves
them through the subsystem's `resolve_edges` helper, which indexes `g`'s edges so
each lookup is a hash probe rather than a scan.

**Complexity / limits.** `O(V + E)` to build the view and edge index, then `O(1)`
per queried edge. The companion `VertexIndex` does the same for vertices. Returns
`NULL` (unevaluated) if any queried edge is absent.

- `Protected`. A non-graph first argument, or a vertex/edge not in the graph, is
  left unevaluated.
- `EdgeIndex` matches `u -> v` against an undirected edge of an undirected
  graph, as Mathematica does; an undirected edge matches either orientation.

**Attributes:** `Protected`.

## References

**See also:** [VertexIndex](../../graphs/VertexIndex/)

- Source: [`src/graph/gops_edit.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_edit.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`EdgeIndex[g, e]` is the inverse of indexing `EdgeList[g]`: it reports where `e`
sits in the edge list, counting from 1. An undirected edge matches in either
orientation, and the sugar `u -> v` / `u <-> v` is accepted alongside explicit
`DirectedEdge` / `UndirectedEdge`.

An edge that is not in the graph leaves the expression unevaluated rather than
returning `0` or a `Missing` wrapper.
