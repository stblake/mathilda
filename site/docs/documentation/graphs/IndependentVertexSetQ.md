# IndependentVertexSetQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IndependentVertexSetQ[g, vlist] gives True if no two vertices of vlist are adjacent in g.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= IndependentVertexSetQ[CycleGraph[6], {1, 3, 5}]
Out[1]= True

In[2]:= IndependentVertexSetQ[CycleGraph[6], {1, 2}]
Out[2]= False

In[3]:= IndependentVertexSetQ[CycleGraph[6], {1, 1, 3}]
Out[3]= True

In[4]:= IndependentVertexSetQ[CycleGraph[6], {1, 7}]
Out[4]= False

In[5]:= VertexCoverQ[CycleGraph[4], {1, 3}]
Out[5]= True

In[6]:= VertexCoverQ[x, {1}]
Out[6]= False
```

### Applications (3)

Opposite corners of a 4-cycle are non-adjacent

```mathematica
In[7]:= IndependentVertexSetQ[CycleGraph[4], {1, 3}]
Out[7]= True
```

Adjacent vertices -- not independent

```mathematica
In[8]:= IndependentVertexSetQ[CycleGraph[4], {1, 2}]
Out[8]= False
```

A single vertex is always independent

```mathematica
In[9]:= IndependentVertexSetQ[CompleteGraph[4], {2}]
Out[9]= True
```

## Implementation notes

**Algorithm.** `IndependentVertexSetQ[g, vlist]` is `True` when no two vertices of
`vlist` are adjacent in `g`. It and `VertexCoverQ` share one routine,
`gm_vertex_set_q` (called with `cover = 0`): `gm_vertex_positions` resolves each
element of `vlist` to a vertex index (returning `False` if any element is not a
vertex of `g`), a membership flag array `in[]` marks the chosen vertices, and a
single pass over the edge-index arrays `(eu, ev)` reports `False` as soon as an
edge has **both** endpoints in the set. Edge direction is ignored (Wolfram's
semantics), so an `UndirectedEdge` and a `DirectedEdge` both count as an
adjacency, and repeated vertices in `vlist` are allowed.

**Data structures.** The edge endpoints come straight from the graph's cached
index arrays via `graph_edge_indices` — no `GalgUG` or bitset is built, since the
predicate only needs one linear edge scan. Membership is a `char in[n]` bitmap
and the chosen positions are a small `int` array; both are freed on every path.

**Complexity / limits.** `O(|vlist| + m)` time and `O(n)` space: one scan to mark
the set, one scan over the edges. The predicate is total — it returns a Boolean
for a non-graph first argument or a non-vertex list element (both `False`) and
never stays unevaluated, as the `*Q` convention requires.

- `Protected` membership predicates; give `False` (never stay unevaluated) when
  `g` is not a graph.
- An element that is not a vertex of `g` gives `False`.
- Repeated vertices are allowed in `vs`.

**Attributes:** `Protected`.

## References

**See also:** [VertexCoverQ](../../graphs/VertexCoverQ/)

- Source: [`src/graph/galg_mis.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_mis.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

`IndependentVertexSetQ[g, vlist]` is `True` when no two vertices of `vlist` are
adjacent in `g`. Edge direction is ignored, so an `UndirectedEdge` and a
`DirectedEdge` both count as an adjacency, and repeated vertices in `vlist` are
allowed.

It is a total predicate: a non-graph first argument, or a `vlist` element that is
not a vertex of `g`, gives `False` rather than staying unevaluated. It is the
complement of `VertexCoverQ` — `vlist` is independent iff its complement is a
vertex cover.
