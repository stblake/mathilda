### Worked examples

```mathematica
In[1]:= EdgeList[AdjacencyGraph[{{0, 1, 0}, {1, 0, 1}, {0, 1, 0}}]]  (* a symmetric matrix gives an undirected graph *)
```

```mathematica
In[1]:= DirectedGraphQ[AdjacencyGraph[{{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}]]  (* an asymmetric matrix gives a directed one *)
```

```mathematica
In[1]:= AdjacencyMatrix[AdjacencyGraph[{{0, 1, 1}, {1, 0, 1}, {1, 1, 0}}]]  (* round-trips with AdjacencyMatrix *)
```

```mathematica
In[1]:= VertexCount[AdjacencyGraph[{{0, 1}, {1, 0}}]]  (* vertices are the integers 1..n *)
```

### Notes

The matrix must be square with every entry a literal `0` or `1`; anything else
leaves the call unevaluated. Diagonal entries are self-loops and are ignored,
so a graph with a nonzero diagonal cannot be represented — Mathilda graphs are
loop-free.

The symmetry of the matrix decides direction: a symmetric matrix becomes an
undirected graph (one edge per unordered pair), an asymmetric one a directed
graph (one arc per nonzero off-diagonal cell). This is the exact inverse of
`AdjacencyMatrix` whenever the source graph's vertices are `1, ..., n`.
