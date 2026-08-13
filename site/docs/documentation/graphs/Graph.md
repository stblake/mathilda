# Graph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Graph[v, e] represents a graph with vertices v and edges e. Graph[e] derives the vertices from the edge list. Edges are DirectedEdge[u,v] or UndirectedEdge[u,v]; u->v and u<->v are accepted as shorthand. Simple graphs only: no self-loops or parallel edges.`**

## Examples

_No verified examples yet for this function._

## Algorithm

construct.c - builtin_graph: normalize, derive, validate, canonicalize.

Accepts:

```text
  Graph[edges]          -- vertices derived from the edges (directed default)
  Graph[verts, edges]   -- explicit vertex list
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

**Attributes:** `Protected`.

## References

**See also:** [Rule](../../assignment-and-rules/Rule/), [InputForm](../../expression-information/InputForm/), [FullForm](../../expression-information/FullForm/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
