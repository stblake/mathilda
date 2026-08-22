# AdjacencyMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AdjacencyMatrix[g] gives the 0/1 adjacency matrix of g (symmetric for undirected graphs).`**

## Examples

_No verified examples yet for this function._

## Algorithm

adjmat.c - AdjacencyMatrix[g]: dense 0/1 adjacency matrix.

Returns an n x n dense List-of-Lists (n = |V|, in canonical vertex order), consumable directly by Det, Tr, and Eigenvalues with no linalg changes. A DirectedEdge[a,b] sets M[a][b] = 1; an UndirectedEdge sets both M[a][b] and M[b][a], so undirected graphs yield a symmetric matrix. Entries are always 0/1 since parallel edges are forbidden.

Future hook: a WeightedAdjacencyMatrix would fill entries with edge weights instead of 1 (Locked Decision 2); not implemented in the MVP.

Memory (SPEC section 4): returns a freshly-allocated matrix; frees res.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
