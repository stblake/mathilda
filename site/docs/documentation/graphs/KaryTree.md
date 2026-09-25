# KaryTree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KaryTree[n] gives a binary tree with n vertices; KaryTree[n, k] a k-ary tree with n vertices, vertices numbered in breadth-first order.`**

## Examples (6)

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

## Implementation notes

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

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
