# AdjacencyMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AdjacencyMatrix[g] gives the 0/1 adjacency matrix of g (symmetric for undirected graphs).`**

## Examples (11)

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

### Applications (6)

Symmetric for an undirected graph

```mathematica
In[6]:= AdjacencyMatrix[CycleGraph[4]]
Out[6]= {{0, 1, 0, 1}, {1, 0, 1, 0}, {0, 1, 0, 1}, {1, 0, 1, 0}}
```

A directed edge sets only M[a][b]

```mathematica
In[7]:= AdjacencyMatrix[Graph[{1 -> 2, 2 -> 3}]]
Out[7]= {{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}
```

Hub row and column

```mathematica
In[8]:= AdjacencyMatrix[StarGraph[4]]
Out[8]= {{0, 1, 1, 1}, {1, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 0, 0}}
```

Six times the number of triangles

```mathematica
In[9]:= Tr[MatrixPower[AdjacencyMatrix[CompleteGraph[3]], 3]]
Out[9]= 6
```

Feeds linear algebra directly

```mathematica
In[10]:= Det[AdjacencyMatrix[CycleGraph[4]]]
Out[10]= 0
```

Rows follow the vertex list order

```mathematica
In[11]:= AdjacencyMatrix[Graph[{a, b, c}, {a <-> c}]]
Out[11]= {{0, 0, 1}, {0, 0, 0}, {1, 0, 0}}
```

## Algorithm

adjmat.c - AdjacencyMatrix[g]: dense 0/1 adjacency matrix.

Returns an n x n dense List-of-Lists (n = |V|, in canonical vertex order), consumable directly by Det, Tr, and Eigenvalues with no linalg changes. A DirectedEdge[a,b] sets M[a][b] = 1; an UndirectedEdge sets both M[a][b] and M[b][a], so undirected graphs yield a symmetric matrix. Entries are always 0/1 since parallel edges are forbidden.

Future hook: a WeightedAdjacencyMatrix would fill entries with edge weights instead of 1 (Locked Decision 2); not implemented in the MVP.

Memory (SPEC section 4): returns a freshly-allocated matrix; frees res.

## Implementation notes

**Algorithm.** `builtin_adjacency_matrix` takes exactly one argument and declines (`NULL`) unless `graph_is_valid` accepts it. It allocates an `n x n` zeroed `int` grid, resolves each edge endpoint to a vertex position, and sets `M[a][b] = 1` for a `DirectedEdge[a, b]`. An `UndirectedEdge` also sets `M[b][a]`, so an undirected graph gives a symmetric matrix. Rows and columns follow the canonical `VertexList` order.

**Data structures.** The graph is `Graph[List[v...], List[edge...]]`. Endpoints resolve through a `GraphVIdx` hash index built once per call. This replaced a linear scan per edge that cost `O(E V)` on top of the matrix itself. The result is a freshly built `List` of `List`s of machine integers, so `Det`, `Tr` and `Eigenvalues` consume it directly. Parallel edges are forbidden by validation, so every entry is `0` or `1`.

**Complexity / limits.** `O(V^2 + E)` time and space, and the `V^2` output dominates. The matrix is dense, not sparse. A non-graph argument, or a call with more or fewer than one argument, is left unevaluated. The result is not cached.

- `Protected`. Symmetric for undirected graphs; entry `(i, j)` is `1` for a
  directed edge `vi -> vj`.
- Linear-algebra interop: the result is an ordinary matrix, so `Det`, `Tr`,
  `Eigenvalues`, `MatrixPower`, etc. apply directly. Inverse: `AdjacencyGraph`.
  Weighted variant: `WeightedAdjacencyMatrix`.
- Mathematica returns a `SparseArray`; Mathilda returns a dense list of lists.
- Unevaluated on a non-graph.

**Attributes:** `Protected`.

## References

**See also:** [Det](../../linear-algebra/Det/), [Tr](../../linear-algebra/Tr/), [Eigenvalues](../../linear-algebra/Eigenvalues/), [MatrixPower](../../linear-algebra/MatrixPower/), [AdjacencyGraph](../../graphs/AdjacencyGraph/), [WeightedAdjacencyMatrix](../../graphs/WeightedAdjacencyMatrix/), [SparseArray](../../data-structures/SparseArray/)

- Source: [`src/graph/adjmat.c`](https://github.com/stblake/mathilda/blob/main/src/graph/adjmat.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

`AdjacencyMatrix` returns a dense `List` of `List`s of `0` and `1`, with rows in `VertexList` order. An undirected edge fills both `M[a][b]` and `M[b][a]`, so such graphs give a symmetric matrix. A directed edge fills only `M[a][b]`. Because the result is an ordinary matrix, `Det`, `Tr`, `MatrixPower` and `Eigenvalues` apply directly.
