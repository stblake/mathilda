# PathGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PathGraphQ[g] gives True if g is a path: connected with at least one vertex and, if undirected, every degree at most 2 and one edge fewer than vertices; if directed, every in- and out-degree at most 1. Mixed graphs are never paths.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= {PathGraphQ[PathGraph[Range[4]]], PathGraphQ[CycleGraph[3]], PathGraphQ[StarGraph[4]], PathGraphQ[Graph[{1->2,2->3}]], PathGraphQ[Graph[{1->2,3->2}]]}
Out[1]= {True, True, False, True, False}

In[2]:= {PathGraphQ[Graph[{1->2,2<->3}]], PathGraphQ[Graph[{1<->2,3<->4}]], PathGraphQ[Graph[{1},{}]], PathGraphQ[Graph[{},{}]], PathGraphQ[x]}
Out[2]= {False, False, True, False, False}
```

## Implementation notes

- `Protected`. A non-graph argument gives `False`.
- Mathematica's definition: at least one vertex, connected, and every degree
  `<= 2` (undirected) or every in/out-degree `<= 1` (directed). So cycles count:
  `PathGraphQ[CycleGraph[3]]` is `True`, as in Mathematica.
- Mixed graphs are never paths.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
