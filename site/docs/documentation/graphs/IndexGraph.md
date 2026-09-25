# IndexGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IndexGraph[g] replaces the vertices of g by 1, 2, ..., n (in VertexList order); IndexGraph[g, r] by r, r+1, ..., r+n-1. Edge weights are kept.`**

## Examples (3)

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

## Implementation notes

- `Protected`. A non-graph argument is left unevaluated.
- The default start is `r = 1`. Edge order and weights are kept.

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
