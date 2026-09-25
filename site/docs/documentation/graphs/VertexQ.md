# VertexQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexQ[g, v] gives True if v is a vertex of the graph g (compared structurally, as by SameQ), and False otherwise.`**

## Examples (4)

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

- `Protected`. Vertices are compared structurally (as by `SameQ`): the vertex
  `2` is not matched by `2.0`. `False` for a non-graph (see
  `UndirectedGraphQ`).

**Attributes:** `Protected`.

## References

**See also:** [SameQ](../../comparisons/SameQ/), [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
