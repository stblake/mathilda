# EdgeRules

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeRules[g] gives the edges of g as a list of rules u -> v (undirected edges too), in EdgeList order.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= EdgeRules[CycleGraph[3]]
Out[1]= {1 -> 2, 2 -> 3, 3 -> 1}

In[2]:= EdgeRules[Graph[{1->2, 2->3}]]
Out[2]= {1 -> 2, 2 -> 3}

In[3]:= EdgeRules[5]
Out[3]= EdgeRules[5]
```

## Implementation notes

- `Protected`. A non-graph argument is left unevaluated.
- Undirected and directed edges both become rules, in `EdgeList` order.

**Attributes:** `Protected`.

## References

**See also:** [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
