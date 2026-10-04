# EdgeDelete

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeDelete[g, e] removes the edge e from g; EdgeDelete[g, {e1, ...}] removes several (each must be an edge of g, else the call stays unevaluated; an undirected edge matches either orientation); EdgeDelete[g, patt] removes every edge matching patt. Vertices, order and the remaining weights are kept.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= EdgeList[EdgeDelete[CycleGraph[4], 1<->2]]
Out[1]= {2 <-> 3, 3 <-> 4, 4 <-> 1}

In[2]:= EdgeList[EdgeDelete[CycleGraph[4], {2<->1, 3<->4}]]
Out[2]= {2 <-> 3, 4 <-> 1}

In[3]:= EdgeList[EdgeDelete[CycleGraph[4], _[1, _]]]
Out[3]= {2 <-> 3, 3 <-> 4, 4 <-> 1}

In[4]:= EdgeDelete[CycleGraph[4], 1->2]
Out[4]= EdgeDelete[Graph[<4 vertices, 4 edges>], 1 -> 2]

In[5]:= EdgeDelete[CycleGraph[4], 1<->3]
Out[5]= EdgeDelete[Graph[<4 vertices, 4 edges>], TwoWayRule[1, 3]]
```

### Applications (1)

Drop the middle edge

```mathematica
In[6]:= EdgeList[EdgeDelete[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4}], 2 <-> 3]]
Out[6]= {1 <-> 2, 3 <-> 4}
```

## Implementation notes

**Algorithm.** `builtin_edge_delete` removes matching edges and rebuilds the graph on the *same*
vertex set (no vertices are dropped). Each argument item is classified individually: an item
containing a pattern is matched against every edge and clears the keep-flag of each hit; a
literal edge is resolved to its `EdgeList` index (and the call stays unevaluated if that edge is
not in the graph). An undirected edge matches either orientation, because edge keys canonicalise
undirected endpoints. The surviving edges and all vertices are copied into a fresh canonical
`Graph`, remapping endpoint indices.

**Data structures.** A `GopsView` exposes the raw vertex/edge arrays plus integer endpoints, so
the literal-edge pass is an integer pass against a `GopsKeySet` (an open-addressing hash of the
argument's edges, probed once per graph edge). The view is taken *after* any pattern matching,
since matching can evaluate arbitrary user code. `EdgeWeight` entries are kept aligned with the
surviving edges.

**Complexity / limits.** `O(V + E)` plus one hash probe per edge for the literal path; the
pattern path costs `O(items·E)` matches. The result is a canonical `Graph`; a non-graph argument,
or a literal that is not an edge of the graph, returns unevaluated.

- `Protected`. A non-graph first argument is left unevaluated.
- Each listed edge must exist, otherwise the call is left unevaluated. An
  undirected edge matches either orientation; `1 -> 2` is not an edge of an
  undirected graph.
- Orders and the remaining weights are kept.
- `O(V + E)` plus one hash per argument item; result seeded into the graph memo.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gops_edit.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_edit.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`EdgeDelete[g, e]` removes the edge `e`, keeping every vertex (the vertex set is unchanged).
`EdgeDelete[g, {e1, ...}]` removes several, and `EdgeDelete[g, patt]` removes every edge matching
a pattern. An undirected edge matches either orientation.

A literal edge that is not in `g` leaves the call unevaluated. The result is a canonical `Graph`
object — inspect it with `EdgeList` or `EdgeCount`.
