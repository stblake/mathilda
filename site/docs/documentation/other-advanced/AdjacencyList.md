# AdjacencyList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AdjacencyList[g] gives the adjacency list of g; AdjacencyList[g,v] gives the vertices adjacent to v (successors for directed edges).`**

## Examples

_No verified examples yet for this function._

## Algorithm

adjlist.c - AdjacencyList[g] and AdjacencyList[g, v].

Neighbors of v: successors for directed edges (v -> u yields u), and both endpoints for undirected edges. This matches the row convention of AdjacencyMatrix (a 1 in row v, column u means v is adjacent to u). Neighbors are returned in first-appearance order, de-duplicated.

```text
AdjacencyList[g]     -> {neighbors(v1), neighbors(v2), ...} in vertex order.
AdjacencyList[g, v]  -> neighbors(v).
```

Memory (SPEC section 4): returns freshly-allocated lists; evaluator frees res.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
