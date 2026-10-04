# VertexDelete

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexDelete[g, v] removes the vertex v and its incident edges from g; VertexDelete[g, {v1, ...}] removes several (each must be a vertex of g, else the call stays unevaluated); VertexDelete[g, patt] removes every vertex matching patt. Vertex and edge order and edge weights are kept.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= EdgeList[VertexDelete[CycleGraph[5], {1, 3}]]
Out[1]= {4 <-> 5}

In[2]:= VertexList[VertexDelete[PathGraph[Range[6]], _?EvenQ]]
Out[2]= {1, 3, 5}

In[3]:= VertexDelete[CycleGraph[3], 7]
Out[3]= VertexDelete[Graph[<3 vertices, 3 edges>], 7]
```

### Options (1)

```mathematica
In[4]:= InputForm[VertexDelete[Graph[{1,2,3,4},{1<->2,2<->3,3<->4},EdgeWeight->{5,6,7}], 2]]
Out[4]= Graph[{1, 3, 4}, {3 <-> 4}, EdgeWeight -> {7}]
```

### Applications (2)

```mathematica
In[5]:= vd = VertexDelete[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}], 3];
```

Vertex 3 and its three incident edges are gone

```mathematica
In[6]:= {VertexList[vd], EdgeList[vd]}
Out[6]= {{1, 2, 4}, {1 -> 2}}
```

## Implementation notes

**Algorithm.** `builtin_vertex_delete` removes one or more vertices together with every incident
edge. It starts with all vertices kept, then clears the keep-flag of each deleted vertex: a
single literal vertex, each item of a list (a literal vertex, else a pattern matched against all
vertices), or a bare pattern matched against all vertices. An edge survives only when *both* of
its endpoints survive. The surviving vertices and edges are copied into a fresh canonical
`Graph`, preserving the original order and remapping endpoint indices.

**Data structures.** A `GopsView` over the raw arrays with integer endpoints, a `vkeep[]` vertex
mask and a derived `ekeep[]` edge mask (`ekeep[k] = vkeep[eu[k]] && vkeep[ev[k]]`). The view is
taken after pattern matching, which can evaluate arbitrary user code.

**Complexity / limits.** `O(V + E)`; the pattern path adds `O(V)` matches per pattern. The result
is a canonical `Graph`. A list item that is neither a vertex of the graph nor a pattern, or a
non-graph argument, returns unevaluated.

- `Protected`. A non-graph first argument is left unevaluated.
- Every listed vertex must exist; otherwise the call is left unevaluated.
- Orders and weights of the surviving vertices and edges are kept.
- `O(V + E)` integer pass over the memoized endpoint arrays plus one hash per
  argument item; the result is seeded into the graph memo (see `VertexAdd`).

**Attributes:** `Protected`.

## References

**See also:** [VertexAdd](../../graphs/VertexAdd/)

- Source: [`src/graph/gops_edit.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_edit.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`VertexDelete[g, v]` removes vertex `v` and every edge touching it. `VertexDelete[g, {v1, ...}]`
removes several, and `VertexDelete[g, patt]` removes every vertex matching a pattern; the
original vertex and edge order is preserved among the survivors.

Each listed vertex must be a vertex of `g` (or the argument must be a pattern), otherwise the call
stays unevaluated. The result is a canonical `Graph` — inspect it with `VertexList` / `EdgeList`.
