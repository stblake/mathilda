# MixedGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MixedGraphQ[g] gives True if g has both directed and undirected edges.`**

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= {MixedGraphQ[Graph[{1->2, 2<->3}]], MixedGraphQ[Graph[{1->2}]], MixedGraphQ[CycleGraph[3]], MixedGraphQ[x]}
Out[1]= {True, False, False, False}
```

## Implementation notes

- `Protected`. A non-graph argument gives `False`.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
