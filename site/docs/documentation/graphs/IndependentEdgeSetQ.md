# IndependentEdgeSetQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IndependentEdgeSetQ[g, elist] gives True if elist is a set of edges of g no two of which share a vertex.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= IndependentEdgeSetQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[4,3]}]
Out[1]= True

In[2]:= IndependentEdgeSetQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[2,3]}]
Out[2]= False

In[3]:= IndependentEdgeSetQ[Graph[{1->2,3->4}], {DirectedEdge[2,1]}]
Out[3]= False

In[4]:= EdgeCoverQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[3,4]}]
Out[4]= True

In[5]:= EdgeCoverQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[2,3]}]
Out[5]= False

In[6]:= EdgeCoverQ[CycleGraph[4], {UndirectedEdge[1,3]}]
Out[6]= False
```

### Applications (2)

A matching: no shared vertex

```mathematica
In[7]:= IndependentEdgeSetQ[CycleGraph[4], {1 <-> 2, 3 <-> 4}]
Out[7]= True
```

Both touch vertex 2

```mathematica
In[8]:= IndependentEdgeSetQ[CycleGraph[4], {1 <-> 2, 2 <-> 3}]
Out[8]= False
```

## Implementation notes

**Algorithm.** `builtin_independent_edge_set_q` tests whether a list of edges is
a **matching** of `g`: a set of edges of `g`, no two sharing a vertex. The shared
`gm_edge_set_q(res, 0)` walks the list, and for each element checks via
`gm_edge_of` that it really is an edge of `g` (accepting `DirectedEdge`/`Rule` and
`UndirectedEdge`/`TwoWayRule`, verified with `graph_has_edge`) and resolves its
two endpoints to vertex indices. A `hit[]` bitmap over the vertices records used
endpoints; if a new edge touches an already-used vertex the answer is `False`. An
element that is not an edge of `g` also gives `False`. (The same routine with the
cover flag set implements `EdgeCoverQ`.)

**Data structures.** A single `char hit[]` array of length `n` (the vertex
count) and the validated-graph membership test `graph_has_edge`; no adjacency is
built. Vertices are resolved with `galg_vertex_arg`.

**Complexity / limits.** `O(k)` edge checks for a `k`-element list, each an `O(1)`
membership probe on the graph memo. A non-graph or non-list argument gives
`False`. `FindIndependentEdgeSet` constructs a *maximum* matching.

- `Protected` membership predicates; give `False` when `g` is not a graph.
- An element that is not an edge of `g` gives `False`. An `UndirectedEdge`
  matches either orientation; a `DirectedEdge` matches only as given.

**Attributes:** `Protected`.

## References

**See also:** [EdgeCoverQ](../../graphs/EdgeCoverQ/)

- Source: [`src/graph/galg_mis.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_mis.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

An independent edge set is a matching: a set of edges of `g` that pairwise share
no vertex. The test verifies every listed element is actually an edge of `g`
(an undirected edge matches either orientation) and that no vertex is used
twice.

A list containing something that is not an edge of `g` gives `False`. Use
`FindIndependentEdgeSet` to construct a maximum matching.
