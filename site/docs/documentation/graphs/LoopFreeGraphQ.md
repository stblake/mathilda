# LoopFreeGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LoopFreeGraphQ[g] gives True if g is a graph with no self-loops -- every valid Mathilda graph -- and False otherwise.`**

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= {SimpleGraphQ[CycleGraph[3]], LoopFreeGraphQ[CycleGraph[3]], SimpleGraphQ[5], LoopFreeGraphQ[{1}]}
Out[1]= {True, True, False, False}
```

## Implementation notes

- `Protected`. A non-graph argument gives `False`.
- Mathilda graphs are always simple, so both give `True` for every valid graph.

**Attributes:** `Protected`.

## References

**See also:** [SimpleGraphQ](../../graphs/SimpleGraphQ/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
