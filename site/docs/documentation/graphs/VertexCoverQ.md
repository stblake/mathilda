# VertexCoverQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexCoverQ[g, vlist] gives True if every edge of g has an endpoint in vlist.`**

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

Two opposite vertices touch all four edges

```mathematica
In[7]:= VertexCoverQ[CycleGraph[4], {1, 3}]
Out[7]= True
```

The opposite edge is left uncovered

```mathematica
In[8]:= VertexCoverQ[CycleGraph[4], {1, 2}]
Out[8]= False
```

Every edge of the triangle meets {1,2}

```mathematica
In[9]:= VertexCoverQ[CompleteGraph[3], {1, 2}]
Out[9]= True
```

## Implementation notes

**Algorithm.** `VertexCoverQ[g, vlist]` is `True` when every edge of `g` has at
least one endpoint in `vlist`. It is the complement test of
`IndependentVertexSetQ` and shares the same routine `gm_vertex_set_q` (called
with `cover = 1`): `gm_vertex_positions` resolves the list to vertex indices
(`False` if any element is not a vertex of `g`), the chosen vertices are marked in
a flag array, and one pass over the edge-index arrays `(eu, ev)` reports `False`
as soon as an edge has **neither** endpoint in the set. Edge direction is ignored
and repeated vertices are allowed, matching Wolfram's semantics. (A vertex cover
and a maximum independent set are exact complements, which is why
`FindVertexCover` in the same file delegates to the independent-set solver and
returns the complement.)

**Data structures.** Endpoints are read from the graph's cached index arrays via
`graph_edge_indices`; membership is a `char in[n]` bitmap and the resolved
positions a small `int` array. No `GalgUG` or bitset is built — a single linear
edge scan suffices — and all scratch is freed on every path.

**Complexity / limits.** `O(|vlist| + m)` time, `O(n)` space. The predicate is
total: it returns a Boolean for a non-graph argument or a non-vertex list element
(both `False`) and never stays unevaluated.

- `Protected` membership predicates; give `False` (never stay unevaluated) when
  `g` is not a graph.
- An element that is not a vertex of `g` gives `False`.
- Repeated vertices are allowed in `vs`.

**Attributes:** `Protected`.

## References

**See also:** [IndependentVertexSetQ](../../graphs/IndependentVertexSetQ/)

- Source: [`src/graph/galg_mis.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_mis.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

`VertexCoverQ[g, vlist]` is `True` when every edge of `g` has at least one
endpoint in `vlist`. Edge direction is ignored and repeated vertices are
allowed.

It is a total predicate — a non-graph argument or a non-vertex element gives
`False`, never an unevaluated result — and it is the exact complement of
`IndependentVertexSetQ`: `vlist` covers `g` iff the vertices outside `vlist` form
an independent set.
