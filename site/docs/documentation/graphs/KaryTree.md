# KaryTree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KaryTree[n] gives a binary tree with n vertices; KaryTree[n, k] a k-ary tree with n vertices, vertices numbered in breadth-first order.`**

## Examples (10)

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

### Applications (4)

A complete binary tree: 1 has children 2,3; 2 has 4,5; 3 has 6,7

```mathematica
In[7]:= EdgeList[KaryTree[7]]
Out[7]= {1 <-> 2, 1 <-> 3, 2 <-> 4, 2 <-> 5, 3 <-> 6, 3 <-> 7}
```

Exactly n vertices

```mathematica
In[8]:= VertexCount[KaryTree[10]]
Out[8]= 10
```

Ternary: vertex 1 has children 2,3,4

```mathematica
In[9]:= EdgeList[KaryTree[7, 3]]
Out[9]= {1 <-> 2, 1 <-> 3, 1 <-> 4, 2 <-> 5, 2 <-> 6, 2 <-> 7}
```

Root and internal nodes have higher degree than the leaves

```mathematica
In[10]:= VertexDegree[KaryTree[7]]
Out[10]= {2, 3, 3, 1, 1, 1, 1}
```

## Implementation notes

**Algorithm.** `builtin_kary_tree` builds a `k`-ary tree on `n` vertices with the
vertices numbered in **breadth-first** order: `KaryTree[n]` is binary (`k = 2`),
`KaryTree[n, k]` is `k`-ary. The helper `kary_tree(n, k)` emits, for each vertex
`i` (0-based internally), the edges to its children `k*i + 1, ..., k*i + k` as
long as the child index is `< n`. This is the implicit-heap layout, so the parent
of vertex `c` is `(c - 1) / k`; the first vertices fill complete levels and the
last level is filled left to right. (The sibling head `CompleteKaryTree[n, k]`
builds the *complete* `k`-ary tree of `n` levels by summing the level sizes and
calling the same helper.)

**Data structures.** A growable `Pairs` buffer of integer endpoint pairs; the
finished graph is assembled by `pairs_graph(n, &p)`, which produces the canonical
`Graph` on vertices `1..n`, undirected, and seeds the memo. No adjacency or
hashing is needed — the tree structure is arithmetic on the indices.

**Complexity / limits.** `O(n)` — one edge per non-root vertex. `n` is bounded by
`GEN_MAX_VERTICES`; a negative `n` or `k`, or an over-cap size, leaves the call
unevaluated.

- Undirected on `1..N`, root `1`, vertices numbered level by level; the edge
  list is the sorted list of pairs `{i, j}`, `i < j`, identical to Wolfram's
  `EdgeList`.
- `KaryTree` counts vertices; `CompleteKaryTree` counts levels.
- Options (`DirectedEdges`, layout options) are not supported.
- Resource limit: more than 10^8 vertices or 5×10^7 edges is left
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [CompleteKaryTree](../../graphs/CompleteKaryTree/), [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gmet_generators.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_generators.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

`KaryTree[n]` is a binary tree on `n` vertices; `KaryTree[n, k]` a `k`-ary tree.
Vertices are numbered breadth-first, so vertex `i` has children `k(i-1)+2, ...,
k(i-1)+k+1` (1-based) up to `n`, and the last level fills left to right.

For the fully-filled tree given a number of *levels* rather than vertices, use
`CompleteKaryTree`.
