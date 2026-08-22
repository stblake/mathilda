# VertexConnectivity

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VertexConnectivity[g] gives the minimum number of vertices whose removal disconnects g.`**

## Examples

_No verified examples yet for this function._

## Algorithm

connectivity.c - ConnectedGraphQ[g] and VertexConnectivity[g].

Both operate on the underlying undirected graph.

```text
  ConnectedGraphQ[g]    True iff g has >= 1 vertex and forms a single
                        connected component.
  VertexConnectivity[g] the least number of vertices whose removal
                        disconnects g (n-1 for a complete graph, 0 if already
                        disconnected or trivial). Computed by brute-force
                        search over vertex subsets -- exact, intended for the
                        small graphs of a pico-CAS.
```

Memory (SPEC section 4): returns freshly-allocated results; frees res.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
