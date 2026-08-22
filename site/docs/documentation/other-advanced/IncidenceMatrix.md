# IncidenceMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IncidenceMatrix[g] gives the vertex-edge incidence matrix of g (oriented: -1 tail, +1 head for directed edges).`**

## Examples

_No verified examples yet for this function._

## Algorithm

incmat.c - IncidenceMatrix[g]: |V| x |E| incidence matrix.

Column j corresponds to edge j (canonical order), row i to vertex i.

```text
  - UndirectedEdge{a,b}: entries (a,j) and (b,j) are 1.
  - DirectedEdge[a,b]:   (a,j) = -1 (tail), (b,j) = 1 (head)  [oriented].
```

Memory (SPEC section 4): returns a freshly-allocated matrix; frees res.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
