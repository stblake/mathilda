# WeaklyConnectedComponents

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`WeaklyConnectedComponents[g] gives the weakly connected components of g.`**

## Examples

_No verified examples yet for this function._

## Algorithm

components.c - connected-component builtins.

```text
  ConnectedComponents[g]          weak components (underlying undirected)
  WeaklyConnectedComponents[g]    same as ConnectedComponents
  StronglyConnectedComponents[g]  strong components (Tarjan) over directed
                                  adjacency; for undirected graphs this
                                  coincides with the weak components.
```

Each returns a List of Lists of vertices, components in first-appearance order, vertices within a component in canonical index order.

Memory (SPEC section 4): returns freshly-allocated lists; frees res.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
