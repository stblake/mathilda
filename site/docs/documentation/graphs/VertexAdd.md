# VertexAdd

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexAdd[g, v] adds the vertex v to the graph g; VertexAdd[g, {v1, v2, ...}] adds several. Vertices already present are ignored; new ones are appended in order. Edge weights are kept.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= VertexList[VertexAdd[CycleGraph[3], 4]]
Out[1]= {1, 2, 3, 4}

In[2]:= VertexList[VertexAdd[CycleGraph[3], {4, 5, 4, 1}]]
Out[2]= {1, 2, 3, 4, 5}

In[3]:= VertexList[VertexAdd[Graph[{a<->b}], {{1,2}}]]
Out[3]= {a, b, {1, 2}}

In[4]:= VertexAdd[x, 4]
Out[4]= VertexAdd[x, 4]
```

### Applications (4)

A single new vertex goes at the end

```mathematica
In[5]:= VertexList[VertexAdd[PathGraph[{1, 2, 3}], 7]]
Out[5]= {1, 2, 3, 7}
```

Existing and repeated vertices are ignored

```mathematica
In[6]:= VertexList[VertexAdd[PathGraph[{1, 2, 3}], {4, 5, 2, 4}]]
Out[6]= {1, 2, 3, 4, 5}
```

New vertices are isolated, the edges are untouched

```mathematica
In[7]:= EdgeCount[VertexAdd[CycleGraph[4], {a, b}]]
Out[7]= 4
```

Works on directed graphs too

```mathematica
In[8]:= VertexCount[VertexAdd[Graph[{1 -> 2}], {x, y, z}]]
Out[8]= 5
```

## Implementation notes

**Algorithm.** `builtin_vertex_add` takes `VertexAdd[g, v]` or `VertexAdd[g, {v1, v2, ...}]`. It first filters the items: any item already a vertex of `g` (a lookup in the memoized vertex index) is dropped, and repeats among the new items are dropped through a scratch hash, so the survivors keep their first-appearance order. If nothing is new, a copy of `g` is returned. Otherwise the new vertices are appended after the existing ones; the edge list, endpoint arrays and any `EdgeWeight` list are carried over unchanged.

**Data structures.** `Graph[List[vertices], List[edges]]` expression tree; the edit reads the memoized endpoint arrays (`gops_view`) and constructs the result with `gops_graph_new`. Vertex and edge nodes are copied by reference count, not deep-copied.

**Complexity / limits.** `O(V + E)` for the copy plus one hash probe per argument item. Added vertices are isolated, so they cannot create loops or parallel edges. A first argument that is not a valid graph, or a call with other than two arguments, is left unevaluated.

- `Protected`. A non-graph first argument is left unevaluated.
- A list always means a list of vertices; to add a vertex that is itself a list,
  wrap it in another list (`VertexAdd[g, {{1, 2}}]`).
- Vertex and edge orders follow Mathematica 15. Existing weights are kept.
- Shared conventions of the graph-editing family (`VertexAdd`, `VertexDelete`,
  `EdgeAdd`, `EdgeDelete`, `Subgraph`, `NeighborhoodGraph`, `VertexReplace`,
  `EdgeRules`, `VertexIndex`, `EdgeIndex`, `IndexGraph`), the transforms, set
  operations, predicates and cycle/path finders: implemented in
  `src/graph/gops_*.c` (header `src/graph/graph_ops.h`). Every head leaves a
  non-graph argument unevaluated (the `*Q` predicates give `False`). Orders —
  of vertices, of edges, of cycle edges — follow Mathematica 15 unless a
  deviation is listed.
- Performance: every edit is an integer pass over the validated-graph memo's
  endpoint arrays: `O(V + E)` plus one hash per argument item, with vertex, edge
  and weight nodes shared into the result. Results are registered with the memo
  (`graph_memo_seed`) from the endpoint arrays already computed and stamped as
  evaluated, so returning a graph costs one vertex hash per vertex — not a full
  re-validation — and the first accessor on the result is a memo hit. At `10^5`
  vertices the edits run in 4–16 ms: 8–200x faster than Mathematica 15 (warm
  and cold) and 20–60x faster than networkx (`benchmarks/93-graph-ops-editing`).
- Deviation (whole editing family): Mathilda graphs are simple, so an edit whose
  result would have a self-loop or parallel edges — `EdgeAdd` of an existing
  edge, `VertexReplace` merging two adjacent vertices — is left unevaluated
  (Mathematica returns a multigraph).

**Attributes:** `Protected`.

## References

**See also:** [VertexDelete](../../graphs/VertexDelete/), [EdgeAdd](../../graphs/EdgeAdd/), [EdgeDelete](../../graphs/EdgeDelete/), [Subgraph](../../graphs/Subgraph/), [NeighborhoodGraph](../../graphs/NeighborhoodGraph/), [VertexReplace](../../graphs/VertexReplace/), [EdgeRules](../../graphs/EdgeRules/), [VertexIndex](../../graphs/VertexIndex/)

- Source: [`src/graph/gops_edit.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_edit.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

Added vertices are appended after the existing ones, in the order given, and never carry edges. A vertex already present is skipped, so `VertexAdd` is idempotent.
