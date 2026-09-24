# WeightedGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`WeightedGraphQ[g] gives True if g carries edge weights (EdgeWeight).`**

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Options (1)

```mathematica
In[1]:= {WeightedGraphQ[Graph[{1,2},{1<->2},EdgeWeight->{3}]], WeightedGraphQ[CycleGraph[3]], EdgeWeightedGraphQ[Graph[{1,2},{1<->2},EdgeWeight->{3}]], EdgeWeightedGraphQ[CycleGraph[3]], WeightedGraphQ[x]}
Out[1]= {True, False, True, False, False}
```

## Implementation notes

- `Protected`. A non-graph argument gives `False`.

**Attributes:** `Protected`.

## References

**See also:** [EdgeWeightedGraphQ](../../graphs/EdgeWeightedGraphQ/), [EdgeWeight](../../graphs/EdgeWeight/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
