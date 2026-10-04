# HypergraphVertexDelete

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphVertexDelete[h, v] or [h, {v1, ...}] deletes the vertices and every hyperedge containing one of them, as VertexDelete does for a Graph.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= InputForm[HypergraphVertexAdd[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {1, 8, 9}]]
Out[1]= Hypergraph[{1, 2, 3, 4, 5, 6, 7, 8, 9}, {{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]

In[2]:= InputForm[HypergraphVertexAdd[Hypergraph[{{1,2}}], {{a,b}}]]
Out[2]= Hypergraph[{1, 2, {a, b}}, {{1, 2}}]

In[3]:= InputForm[HypergraphVertexDelete[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 3]]
Out[3]= Hypergraph[{1, 2, 4, 5, 6, 7}, {{4, 5, 6}, {7}}]

In[4]:= InputForm[HypergraphVertexDelete[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], {1, 7}]]
Out[4]= Hypergraph[{2, 3, 4, 5, 6}, {{3, 4}, {4, 5, 6}}]

In[5]:= HypergraphVertexDelete[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}], 99]
Out[5]= HypergraphVertexDelete[Hypergraph[<7 vertices, 4 hyperedges>], 99]
```

### Applications (3)

```mathematica
In[6]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
Out[6]= Hypergraph[<7 vertices, 4 hyperedges>]
```

Also removes every hyperedge containing 3

```mathematica
In[7]:= InputForm[HypergraphVertexDelete[h, 3]]
Out[7]= Hypergraph[{1, 2, 4, 5, 6, 7}, {{4, 5, 6}, {7}}]
```

Delete several vertices at once

```mathematica
In[8]:= InputForm[HypergraphVertexDelete[h, {1, 7}]]
Out[8]= Hypergraph[{2, 3, 4, 5, 6}, {{3, 4}, {4, 5, 6}}]
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_vertex_delete` removes the named vertices **and
every hyperedge incident to one of them**, exactly as `VertexDelete` does for a
`Graph`. `edit_items` reads the second argument as one vertex or a List of
vertices; each is resolved through the memo's vertex index, and the call is left
unevaluated if a named vertex is absent. A keep-mask `kv` over vertices and a
keep-mask `ke` over hyperedges are formed (a hyperedge is dropped as soon as one
of its distinct members is deleted), and `rebuild` emits the surviving vertices
and hyperedges in their original order.

**Data structures.** The memo's vertex `GraphVIdx`, the distinct-vertex CSR
`soff/sv` (incidence test per hyperedge), and two `unsigned char` keep-masks;
`rebuild` copies the survivors into a fresh `Hypergraph`.

**Complexity / limits.** `O(n + Σ|e|)`. Mutators unshare the memoized node, so the
input object is never changed.

- In the Vertex heads any List is a list of vertices; wrap a List-valued vertex
  as `{{...}}`.
- `HypergraphVertexDelete` is unevaluated if a named vertex is absent.
- Mutators unshare a memoized hypergraph node before editing, so the original
  object is never changed.

**Attributes:** `Protected`.

## References

**See also:** [HypergraphVertexAdd](../../hypergraphs/HypergraphVertexAdd/), [VertexDelete](../../graphs/VertexDelete/)

- Source: [`src/graph/hyp_ops.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_ops.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`HypergraphVertexDelete[h, v]` removes the vertex `v` **and every hyperedge
containing it**, exactly as `VertexDelete` does for a `Graph`; the list form
removes several. Deleting vertex `3` from the running example therefore drops both
hyperedges `{1, 2, 3}` and `{3, 4}`.

A named vertex that is not in `h` leaves the call unevaluated. As with
`HypergraphVertexAdd`, any `List` argument is a list of vertices, and the edit
produces a fresh object without mutating the original. To delete hyperedges while
keeping their vertices, use `HypergraphEdgeDelete`; to intersect hyperedges with a
vertex set instead of dropping them, use `HypergraphRestriction`.
