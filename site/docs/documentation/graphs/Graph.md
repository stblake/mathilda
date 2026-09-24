# Graph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Graph[v, e] represents a graph with vertices v and edges e. Graph[e] derives the vertices from the edge list. Edges are DirectedEdge[u,v] or UndirectedEdge[u,v]; u->v and u<->v are accepted as shorthand. Simple graphs only: no self-loops or parallel edges.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (3)

```mathematica
In[1]:= Graph[{1,2,3,4}, {1->2, 2->3, 3->4, 4->1}]
Out[1]= Graph[<4 vertices, 4 edges>]

In[2]:= InputForm[Graph[{1,2}, {1<->2}]]
Out[2]= Graph[{1, 2}, {1 <-> 2}]

In[3]:= InputForm[Graph[{1->2, 2->3, 3->1}]]
Out[3]= Graph[{1, 2, 3}, {1 -> 2, 2 -> 3, 3 -> 1}]
```

### Scope (7)

```mathematica
In[4]:= InputForm[Graph[{a,b,c}, {DirectedEdge[a,b], UndirectedEdge[b,c]}]]
Out[4]= Graph[{a, b, c}, {a -> b, b <-> c}]

In[5]:= InputForm[Graph[{1,2,3}, {1->2, 2->3}, EdgeWeight -> {5, 7}]]
Out[5]= Graph[{1, 2, 3}, {1 -> 2, 2 -> 3}, EdgeWeight -> {5, 7}]

In[6]:= FullForm[Graph[{1,2},{1<->2}]]
Out[6]= Graph[List[1, 2], List[UndirectedEdge[1, 2]]]

In[7]:= Graph[{1,2}, {1->1}]
Out[7]= Graph[{1, 2}, {1 -> 1}]

In[8]:= Graph[{1,2}, {1->2, 1->2}]
Out[8]= Graph[{1, 2}, {1 -> 2, 1 -> 2}]

In[9]:= Graph[{1,2}, {1->3}]
Out[9]= Graph[{1, 2}, {1 -> 3}]

In[10]:= Graph[{1,2,3}, {1->2, 2->3}, EdgeWeight -> {5}]
Out[10]= Graph[{1, 2, 3}, {1 -> 2, 2 -> 3}, EdgeWeight -> {5}]
```

## Options & behaviour

### Scope

**Malformed input** is returned unevaluated (a self-loop, a duplicate edge, an
endpoint missing from the vertex list, a weight list of the wrong length):

## Algorithm

construct.c - builtin_graph: normalize, derive, validate, canonicalize.

Accepts:

```text
  Graph[edges]                        -- vertices derived from the edges (directed default)
  Graph[verts, edges]                 -- explicit vertex list
  Graph[verts, edges, EdgeWeight -> {w1, ..., wm}]
                                       -- explicit vertex list + per-edge weights, matched
                                          to `edges` by position; wrong length is malformed
                                          (left unevaluated), same as any other rejection
                                          below. Weighted graphs require the explicit-vertex
                                          form -- Graph[edges, EdgeWeight -> {...}] is not
                                          accepted (deliberately out of scope; see the plan).
```

Edge sugar is normalized on construction:

```text
  Rule[u,v]        / u -> v    ->  DirectedEdge[u, v]
  TwoWayRule[u,v]  / u <-> v   ->  UndirectedEdge[u, v]
  DirectedEdge[u,v] / UndirectedEdge[u,v]   pass through unchanged
```

The result is the canonical Graph[List[verts], List[edges]] with vertices in first-appearance order (when derived). Malformed input -- 3-arg edges, self-loops, parallel edges, or an edge endpoint absent from an explicit vertex list -- leaves Graph[...] unevaluated (returns NULL).

Memory (SPEC section 4): the canonical tree is built entirely from expr_copy of the argument's parts, so `res` is never cannibalized. On success the evaluator frees `res`; on NULL it retains it. The "already canonical" case returns NULL so evaluation reaches a fixed point.

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Graph construction, 40000 edges | 5.81 s | 3.98 s | 14.7 s |
| VertexDegree, 20000 vertices | 0 s | 0.045 s | 0.981 s |
| EdgeCount | 0 s | 0 s | 1.07 s |
| ConnectedComponents | 0 s | 0.001 s | 6.39 s |
| GraphDistance from vertex 1 | 0 s | 1.05e+04 s | 3.91 s |
| FindShortestPath 1 to 10000 | 0 s | 0.692 s | 0.057 s |

## Implementation notes

- `Protected`. A graph is a value: the constructor normalizes and validates its
  input and returns the canonical `Graph[List[verts], List[edges]]` (or, when
  weighted, `Graph[List[verts], List[edges], EdgeWeight -> List[weights]]`).
- Edge normalization: `u -> v` (`Rule`) and `DirectedEdge[u, v]` become
  `DirectedEdge[u, v]`; `u <-> v` (`TwoWayRule`) and `UndirectedEdge[u, v]`
  become `UndirectedEdge[u, v]`. Directed and undirected edges may be mixed.
- Malformed input is left unevaluated: self-loops, parallel/duplicate edges,
  3-argument edges, an edge endpoint absent from an explicit vertex list, or
  (for a weighted graph) an `EdgeWeight` list whose length doesn't match the
  edge list. Anti-parallel directed edges `u -> v` and `v -> u` are distinct and
  allowed.
- The weighted form requires the explicit-vertex form;
  `Graph[e, EdgeWeight -> {...}]` (derived vertices) is not accepted and stays
  unevaluated. A weight list whose length doesn't match `e` is malformed, like
  any other rejection above. Read the weights back with `EdgeWeight`.
- Printing: in standard output a graph shows a terse summary,
  `Graph[<n vertices, m edges>]`. `InputForm` and `FullForm` print the literal
  constructor, which round-trips through the parser.
- Validation is memoized per graph node (see the *Performance model* section of
  this file's preamble), so repeated queries on the same graph do not re-check
  it. `GraphQ` tests validity.

**Attributes:** `Protected`.

## References

**See also:** [Rule](../../assignment-and-rules/Rule/), [EdgeWeight](../../graphs/EdgeWeight/), [InputForm](../../expression-information/InputForm/), [FullForm](../../expression-information/FullForm/), [GraphQ](../../graphs/GraphQ/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
