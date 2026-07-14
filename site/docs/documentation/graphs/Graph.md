# Graph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Graph[v, e] represents a graph with vertices v and edges e. Graph[e] derives the vertices from the edge list. Edges are DirectedEdge[u,v] or UndirectedEdge[u,v]; u->v and u<->v are accepted as shorthand. Simple graphs only: no self-loops or parallel edges.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
