# VertexQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexQ[g, v] gives True if v is a vertex of the graph g (compared structurally, as by SameQ), and False otherwise.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= VertexQ[CycleGraph[3], 2]
Out[1]= True

In[2]:= VertexQ[CycleGraph[3], 2.0]
Out[2]= False

In[3]:= VertexQ[CycleGraph[3], 7]
Out[3]= False

In[4]:= VertexQ[5, 1]
Out[4]= False
```

### Applications (4)

2 is one of the vertices 1..4

```mathematica
In[5]:= VertexQ[CompleteGraph[4], 2]
Out[5]= True
```

9 is not

```mathematica
In[6]:= VertexQ[CompleteGraph[4], 9]
Out[6]= False
```

Structural equality: 1.0 is not the vertex 1

```mathematica
In[7]:= VertexQ[CycleGraph[3], 1.0]
Out[7]= False
```

A non-graph first argument gives False, not unevaluated

```mathematica
In[8]:= VertexQ[5, 2]
Out[8]= False
```

## Algorithm

membership.c - VertexQ[g, v] and EdgeQ[g, e].

```text
  VertexQ[g, v]   True iff v is a vertex of g
  EdgeQ[g, e]     True iff e is an edge of g
```

Membership is structural (SameQ, via expr_eq), as in the Wolfram Language: a vertex 1 is not matched by 1.0.

EdgeQ accepts the same edge sugar the constructor does -- u -> v means DirectedEdge[u, v], u <-> v means UndirectedEdge[u, v]. An undirected query matches an UndirectedEdge in either orientation; a directed query matches only a DirectedEdge with the same ordered endpoints. Direction is never blurred: u -> v is not an edge of Graph[{u, v}, {u <-> v}].

Both give False for anything that is not a valid graph. Both are O(1) hash probes into the validated-graph memo (graph_util.c), which already holds the vertex index and edge-key set that validating g built.

Memory (SPEC section 4): returns a fresh symbol; the evaluator frees res.

## Implementation notes

**Algorithm.** `builtin_vertex_q` answers `VertexQ[g, v]` — `True` iff `v` is a
vertex of `g`. It calls `graph_vertex_position(g, v)` and returns the boolean
`position >= 0`. Membership is structural, matching the Wolfram Language's
`SameQ` semantics via `expr_eq`: the integer vertex `1` is *not* matched by the
real `1.0`. A `g` that is not a valid graph makes `graph_vertex_position` return
`-2`, so `VertexQ` gives `False` rather than leaving itself unevaluated. The head
requires exactly two arguments; any other arity returns `NULL`.

**Data structures.** The position lookup is an `O(1)` probe into the
validated-graph memo (`graph_util.c`), which already holds the `GraphVIdx`
`expr_hash` index of vertices that validating `g` built — no per-call scan of the
vertex list. The only allocation is the returned `True`/`False` symbol.

**Complexity / limits.** `O(1)` expected on a memo hit; the memo itself is built
in `O(V + E)` the first time a graph is validated. Because equality is
structural, callers that want numeric-insensitive membership must normalise
their vertices first.

- `Protected`. Vertices are compared structurally (as by `SameQ`): the vertex
  `2` is not matched by `2.0`. `False` for a non-graph (see
  `UndirectedGraphQ`).

**Attributes:** `Protected`.

## References

**See also:** [SameQ](../../comparisons/SameQ/), [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- Source: [`src/graph/membership.c`](https://github.com/stblake/mathilda/blob/main/src/graph/membership.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`VertexQ[g, v]` tests membership structurally, as `SameQ` does: the integer
vertex `1` is not matched by the real `1.0`. The test is an `O(1)` hash probe
into the graph's validated memo.

Unlike most graph accessors, `VertexQ` returns `False` (rather than staying
unevaluated) when the first argument is not a valid graph.
