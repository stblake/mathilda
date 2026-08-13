# VertexInDegree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexInDegree[g] / VertexInDegree[g,v] gives in-degrees (incoming directed edges; undirected edges count for both).`**

## Examples

_No verified examples yet for this function._

## Algorithm

degree.c - VertexDegree, VertexInDegree, VertexOutDegree.

Each accepts VertexDegree[g] (a list of degrees, one per vertex in canonical order) or VertexDegree[g, v] (the degree of a single vertex).

Conventions (documented in docs/spec/builtins/graphs.md):

```text
  - A DirectedEdge[a,b] adds 1 to out(a) and 1 to in(b); total degree of a
    vertex is in + out, so for a purely directed graph total = in + out.
  - An UndirectedEdge[a,b] is incident to both a and b, adding 1 to each of
    their in-, out-, and total degrees (so in = out = total for a purely
    undirected graph).
```

There are no self-loops, so no endpoint is double-counted within one edge.

Memory (SPEC section 4): returns freshly-allocated integers/lists; the evaluator frees res.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
