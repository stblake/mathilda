# CompleteKaryTree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CompleteKaryTree[n] gives the complete binary tree with n levels; CompleteKaryTree[n, k] the complete k-ary tree with n levels.`**

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= EdgeList[KaryTree[7]]
Out[1]= {1 <-> 2, 1 <-> 3, 2 <-> 4, 2 <-> 5, 3 <-> 6, 3 <-> 7}

In[2]:= EdgeList[KaryTree[5, 3]]
Out[2]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 2 <-> 5}

In[3]:= CompleteKaryTree[3]
Out[3]= Graph[<7 vertices, 6 edges>]

In[4]:= CompleteKaryTree[3, 3]
Out[4]= Graph[<13 vertices, 12 edges>]

In[5]:= EdgeList[CompleteKaryTree[2, 3]]
Out[5]= {1 <-> 2, 1 <-> 3, 1 <-> 4}

In[6]:= KaryTree[0]
Out[6]= KaryTree[0]
```

### Applications (6)

Three levels of a binary tree, 7 vertices

```mathematica
In[7]:= EdgeList[CompleteKaryTree[3]]
Out[7]= {1 <-> 2, 1 <-> 3, 2 <-> 4, 2 <-> 5, 3 <-> 6, 3 <-> 7}
```

A root with three children

```mathematica
In[8]:= EdgeList[CompleteKaryTree[2, 3]]
Out[8]= {1 <-> 2, 1 <-> 3, 1 <-> 4}
```

1 + 3 + 9 + 27 vertices

```mathematica
In[9]:= VertexCount[CompleteKaryTree[4, 3]]
Out[9]= 40
```

Root degree 2, internal nodes 3, leaves 1

```mathematica
In[10]:= VertexDegree[CompleteKaryTree[3, 2]]
Out[10]= {2, 3, 3, 1, 1, 1, 1}
```

A tree has one fewer edge than vertices

```mathematica
In[11]:= EdgeCount[CompleteKaryTree[5, 2]]
Out[11]= 30
```

The result is an ordinary Graph expression

```mathematica
In[12]:= GraphQ[CompleteKaryTree[3, 4]]
Out[12]= True
```

## Implementation notes

**Algorithm.** `builtin_complete_kary_tree` accepts `CompleteKaryTree[n]` (arity 2) or `CompleteKaryTree[n, k]`, where `n` is the number of **levels** and `k` the number of children per node. It sums the level sizes `1 + k + k^2 + ...` to get the vertex count, then hands off to `kary_tree`, which numbers vertices in heap order: vertex `i` (0-based) has children `k*i + 1 .. k*i + k`, those below the vertex count being kept.

**Data structures.** Edges go into the shared packed-64-bit `Pairs` buffer and are produced already in sorted order, so `pairs_graph` skips the sort and builds `Graph[Range[N], {UndirectedEdge[parent, child], ...}]` directly. The result is an ordinary `Expr` tree with vertices `1..N`.

**Complexity / limits.** `O(N)` for `N` vertices, with `N - 1` edges. The level-size loop checks overflow and refuses a tree exceeding 10^8 vertices; non-integer, zero or negative arguments leave the call unevaluated. `KaryTree[n, k]` shares the same builder but takes a vertex count instead of a level count.

- Undirected on `1..N`, root `1`, vertices numbered level by level; the edge
  list is the sorted list of pairs `{i, j}`, `i < j`, identical to Wolfram's
  `EdgeList`.
- `KaryTree` counts vertices; `CompleteKaryTree` counts levels.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [KaryTree](../../graphs/KaryTree/), [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gmet_generators.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_generators.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
- Tests: [`tests/test_graphplot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphplot.c)

## Notes & additional examples

### Notes

The first argument counts **levels**, not vertices, so `CompleteKaryTree[n, k]` has `(k^n - 1)/(k - 1)` vertices. The second argument is the number of children per node and defaults to 2. Vertices are numbered in heap order: the children of vertex `i` are `k(i - 1) + 2 .. k(i - 1) + k + 1`.

Use `KaryTree[n, k]` instead when you want a given number of vertices with the last level partly filled.
