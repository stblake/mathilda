# AdjacencyMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AdjacencyMatrix[g] gives the 0/1 adjacency matrix of g (symmetric for undirected graphs).`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= AdjacencyMatrix[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]
Out[1]= {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}

In[2]:= Det[AdjacencyMatrix[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]]
Out[2]= -1

In[3]:= AdjacencyMatrix[PathGraph[3]]
Out[3]= {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}}

In[4]:= Eigenvalues[AdjacencyMatrix[CycleGraph[4]]]
Out[4]= {-2, 2, 0, 0}

In[5]:= AdjacencyMatrix[Graph[{},{}]]
Out[5]= {}
```

## Algorithm

adjmat.c - AdjacencyMatrix[g]: dense 0/1 adjacency matrix.

Returns an n x n dense List-of-Lists (n = |V|, in canonical vertex order), consumable directly by Det, Tr, and Eigenvalues with no linalg changes. A DirectedEdge[a,b] sets M[a][b] = 1; an UndirectedEdge sets both M[a][b] and M[b][a], so undirected graphs yield a symmetric matrix. Entries are always 0/1 since parallel edges are forbidden.

Future hook: a WeightedAdjacencyMatrix would fill entries with edge weights instead of 1 (Locked Decision 2); not implemented in the MVP.

Memory (SPEC section 4): returns a freshly-allocated matrix; frees res.

## Implementation notes

- `Protected`. Symmetric for undirected graphs; entry `(i, j)` is `1` for a
  directed edge `vi -> vj`.
- Linear-algebra interop: the result is an ordinary matrix, so `Det`, `Tr`,
  `Eigenvalues`, `MatrixPower`, etc. apply directly. Inverse: `AdjacencyGraph`.
  Weighted variant: `WeightedAdjacencyMatrix`.
- Mathematica returns a `SparseArray`; Mathilda returns a dense list of lists.
- Unevaluated on a non-graph.

**Attributes:** `Protected`.

## References

**See also:** [Det](../../linear-algebra/Det/), [Tr](../../linear-algebra/Tr/), [Eigenvalues](../../linear-algebra/Eigenvalues/), [MatrixPower](../../linear-algebra/MatrixPower/), [AdjacencyGraph](../../graphs/AdjacencyGraph/), [WeightedAdjacencyMatrix](../../graphs/WeightedAdjacencyMatrix/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
