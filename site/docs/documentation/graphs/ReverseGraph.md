# ReverseGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ReverseGraph[g] reverses every directed edge of g; undirected edges, the edge order and the weights are kept.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= EdgeList[ReverseGraph[Graph[{1->2,2->3,3<->4}]]]
Out[1]= {2 -> 1, 3 -> 2, 3 <-> 4}

In[2]:= EdgeList[ReverseGraph[CycleGraph[3]]]
Out[2]= {1 <-> 2, 2 <-> 3, 3 <-> 1}
```

### Options (1)

```mathematica
In[3]:= InputForm[ReverseGraph[Graph[{1,2,3},{1->2,3->2},EdgeWeight->{4,5}]]]
Out[3]= Graph[{1, 2, 3}, {2 -> 1, 2 -> 3}, EdgeWeight -> {4, 5}]
```

## Implementation notes

- `Protected`. A non-graph argument is left unevaluated.
- Undirected edges are kept as they are; edge order and weights are kept.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
