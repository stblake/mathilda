# PathGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PathGraphQ[g] gives True if g is a path: connected with at least one vertex and, if undirected, every degree at most 2 and one edge fewer than vertices; if directed, every in- and out-degree at most 1. Mixed graphs are never paths.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= {PathGraphQ[PathGraph[Range[4]]], PathGraphQ[CycleGraph[3]], PathGraphQ[StarGraph[4]], PathGraphQ[Graph[{1->2,2->3}]], PathGraphQ[Graph[{1->2,3->2}]]}
Out[1]= {True, True, False, True, False}

In[2]:= {PathGraphQ[Graph[{1->2,2<->3}]], PathGraphQ[Graph[{1<->2,3<->4}]], PathGraphQ[Graph[{1},{}]], PathGraphQ[Graph[{},{}]], PathGraphQ[x]}
Out[2]= {False, False, True, False, False}
```

### Applications (4)

A simple path

```mathematica
In[3]:= PathGraphQ[PathGraph[{1, 2, 3, 4}]]
Out[3]= True
```

A cycle counts as a closed path

```mathematica
In[4]:= PathGraphQ[CycleGraph[3]]
Out[4]= True
```

The hub has degree 3, so not a path

```mathematica
In[5]:= PathGraphQ[StarGraph[4]]
Out[5]= False
```

Too many edges for every degree <= 2

```mathematica
In[6]:= PathGraphQ[CompleteGraph[4]]
Out[6]= False
```

## Implementation notes

**Algorithm.** `builtin_path_graph_q` follows Mathematica's definition of a path:
at least one vertex, connected, and with bounded degrees — undirected, every
degree `<= 2`; directed, every in- and out-degree `<= 1`. One pass over the edges
accumulates in/out degrees and fails as soon as a bound is exceeded; a union-find
over the edges then checks connectivity (the number of joins must equal `n - 1`).
Because a cycle meets the degree bounds and is connected, it counts as a (closed)
path — `PathGraphQ[CycleGraph[3]]` is `True`, as in Mathematica 15. A mixed graph
is never a path, and a graph with more than `n` edges fails the degree bounds
immediately.

**Data structures.** A `GopsView` plus three `O(V)` arrays: in-degree `din[]`,
out-degree `dout[]`, and the union-find `parent[]`. All work is on the view's
integer endpoint arrays.

**Complexity / limits.** `O(V + E)`. A non-graph argument gives `False`. The
early `ne > n` rejection keeps a dense graph from even allocating the degree
arrays.

- `Protected`. A non-graph argument gives `False`.
- Mathematica's definition: at least one vertex, connected, and every degree
  `<= 2` (undirected) or every in/out-degree `<= 1` (directed). So cycles count:
  `PathGraphQ[CycleGraph[3]]` is `True`, as in Mathematica.
- Mixed graphs are never paths.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gops_preds.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_preds.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)

## Notes & additional examples

### Notes

`PathGraphQ[g]` is `True` when `g` is connected and every vertex has degree at
most 2 (in- and out-degree at most 1 for a directed graph). Following
Mathematica, a cycle satisfies this and so is a path; a vertex of degree 3 (as
in a star) or the density of a complete graph on 4+ vertices fails it.

A mixed graph is never a path.
