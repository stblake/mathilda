# AdjacencyGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AdjacencyGraph[m] builds a graph on vertices 1..n from a 0/1 adjacency matrix m (undirected if m is symmetric, else directed).`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= AdjacencyGraph[{{0,1},{1,0}}]
Out[1]= Graph[<2 vertices, 1 edge>]

In[2]:= EdgeList[AdjacencyGraph[{{0,1,0},{0,0,1},{1,0,0}}]]
Out[2]= {1 -> 2, 2 -> 3, 3 -> 1}

In[3]:= EdgeList[AdjacencyGraph[{{0,1,1},{1,0,0},{1,0,0}}]]
Out[3]= {1 <-> 2, 1 <-> 3}

In[4]:= AdjacencyMatrix[AdjacencyGraph[AdjacencyMatrix[CycleGraph[4]]]] == AdjacencyMatrix[CycleGraph[4]]
Out[4]= True

In[5]:= AdjacencyGraph[{{1,2},{3}}]
Out[5]= AdjacencyGraph[{{1, 2}, {3}}]
```

### Applications (4)

A symmetric matrix gives an undirected graph

```mathematica
In[6]:= EdgeList[AdjacencyGraph[{{0, 1, 0}, {1, 0, 1}, {0, 1, 0}}]]
Out[6]= {1 <-> 2, 2 <-> 3}
```

An asymmetric matrix gives a directed one

```mathematica
In[7]:= DirectedGraphQ[AdjacencyGraph[{{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}]]
Out[7]= True
```

Round-trips with AdjacencyMatrix

```mathematica
In[8]:= AdjacencyMatrix[AdjacencyGraph[{{0, 1, 1}, {1, 0, 1}, {1, 1, 0}}]]
Out[8]= {{0, 1, 1}, {1, 0, 1}, {1, 1, 0}}
```

Vertices are the integers 1..n

```mathematica
In[9]:= VertexCount[AdjacencyGraph[{{0, 1}, {1, 0}}]]
Out[9]= 2
```

## Algorithm

adjgraph.c - AdjacencyGraph[m]: build a graph from a 0/1 adjacency matrix.

The inverse of AdjacencyMatrix. Vertices are the integers 1..n. A symmetric matrix yields an undirected graph (one UndirectedEdge per i<j with m[i][j]=1); an asymmetric matrix yields a directed graph (a DirectedEdge for each off- diagonal m[i][j]=1). Diagonal entries (self-loops) are ignored. The result is returned as a Graph[...] expression and canonicalized/validated by the evaluator (builtin_graph).

Round-trips with AdjacencyMatrix when the source graph's vertices are 1..n.

Memory (SPEC section 4): returns a freshly-allocated Graph; frees res.

## Implementation notes

**Algorithm.** `builtin_adjacency_graph` is the inverse of `AdjacencyMatrix`: it
reads an `n x n` matrix of literal `0`/`1` integers and builds a graph on the
vertices `1..n`. It first tests the matrix for symmetry; a symmetric matrix
yields an **undirected** graph with one `UndirectedEdge[i, j]` per pair `i < j`
with `m[i][j] = 1`, and an asymmetric matrix yields a **directed** graph with a
`DirectedEdge[i, j]` for each off-diagonal `m[i][j] = 1`. Diagonal entries
(self-loops) are ignored. Any entry that is not `0` or `1`, or a non-square or
non-list argument, leaves the call unevaluated (returns `NULL`).

**Data structures.** The matrix rows are read directly as `Expr` lists. The
builder allocates the vertex list `{1, ..., n}` and an upper-bounded edge array
(`n*n` slots, trimmed by the actual count), then wraps them in a
`Graph[List, List]` expression. That raw form is handed back to the evaluator,
whose `builtin_graph` normalizes and validates it and seeds the per-node graph
memo (the shared `GraphVIdx`/endpoint arrays described in `src/graph/graph.h`).

**Complexity / limits.** `O(n^2)` — the symmetry test and the fill are both full
passes over the matrix. The round-trip `AdjacencyGraph[AdjacencyMatrix[g]]`
reproduces `g` exactly when `g`'s vertices are `1..n`. Weights are not read: the
input is a plain 0/1 adjacency matrix, so the result is unweighted.

- `Protected`. The inverse of `AdjacencyMatrix`: undirected if `m` is symmetric,
  else directed. `AdjacencyGraph[AdjacencyMatrix[g]]` reproduces `g`'s edges
  (as a set — edges are regenerated in row-major order, so `EdgeList` order and
  undirected-edge orientation may differ from `g`'s; the adjacency matrices
  agree).
- A non-square (ragged) matrix, or entries other than 0/1, leave the call
  unevaluated. A `1` on the diagonal is currently dropped silently rather than
  rejected (the graph model has no self-loops).

**Attributes:** `Protected`.

## References

**See also:** [AdjacencyMatrix](../../graphs/AdjacencyMatrix/), [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/adjgraph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/adjgraph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

The matrix must be square with every entry a literal `0` or `1`; anything else
leaves the call unevaluated. Diagonal entries are self-loops and are ignored,
so a graph with a nonzero diagonal cannot be represented — Mathilda graphs are
loop-free.

The symmetry of the matrix decides direction: a symmetric matrix becomes an
undirected graph (one edge per unordered pair), an asymmetric one a directed
graph (one arc per nonzero off-diagonal cell). This is the exact inverse of
`AdjacencyMatrix` whenever the source graph's vertices are `1, ..., n`.
