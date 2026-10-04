# IndexGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IndexGraph[g] replaces the vertices of g by 1, 2, ..., n (in VertexList order); IndexGraph[g, r] by r, r+1, ..., r+n-1. Edge weights are kept.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= EdgeList[IndexGraph[Graph[{a<->b, b<->c}]]]
Out[1]= {1 <-> 2, 2 <-> 3}

In[2]:= EdgeList[IndexGraph[Graph[{a<->b, b<->c}], 10]]
Out[2]= {10 <-> 11, 11 <-> 12}

In[3]:= IndexGraph[{1,2}]
Out[3]= IndexGraph[{1, 2}]
```

### Applications (4)

Vertices become 1, 2, 3

```mathematica
In[4]:= VertexList[IndexGraph[Graph[{a, b, c}, {a <-> b, b <-> c}]]]
Out[4]= {1, 2, 3}
```

An offset of 0 starts at 0

```mathematica
In[5]:= EdgeList[IndexGraph[Graph[{a, b, c}, {a <-> b, b <-> c}], 0]]
Out[5]= {0 <-> 1, 1 <-> 2}
```

Directed edges keep their direction

```mathematica
In[6]:= EdgeList[IndexGraph[Graph[{x -> y, y -> z}]]]
Out[6]= {1 -> 2, 2 -> 3}
```

Offset labels r, r+1, ...

```mathematica
In[7]:= EdgeList[IndexGraph[CycleGraph[4], 10]]
Out[7]= {10 <-> 11, 11 <-> 12, 12 <-> 13, 13 <-> 10}
```

## Implementation notes

**Algorithm.** `builtin_index_graph` handles `IndexGraph[g]` and `IndexGraph[g, r]` with `r` a machine integer (default 1). Vertex `i` of `VertexList[g]` is renamed to the integer `r + i`, so the vertices become `r, r+1, ...` in their existing order. Every edge is rebuilt from the same endpoint indices with its own directedness, so orientation and edge order are preserved, and any `EdgeWeight` list is copied across. An offset that would overflow `int64` leaves the call unevaluated.

**Data structures.** `Graph[List, List]` expression tree; because the renaming is positional, the endpoint arrays `eu`/`ev`/`directed` of the memo are reused unchanged, and only the vertex and edge expressions are new.

**Complexity / limits.** `O(V + E)`, with no hashing needed. A non-integer `r` or an invalid graph is left unevaluated.

- `Protected`. A non-graph argument is left unevaluated.
- The default start is `r = 1`. Edge order and weights are kept.

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/gops_edit.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_edit.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

Vertex number `i` in `VertexList` is renamed `r + i - 1`, with `r = 1` by default, so the structure of the graph is unchanged and only the labels differ. The offset must be a machine integer.
