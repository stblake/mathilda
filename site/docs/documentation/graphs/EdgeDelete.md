# EdgeDelete

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeDelete[g, e] removes the edge e from g; EdgeDelete[g, {e1, ...}] removes several (each must be an edge of g, else the call stays unevaluated; an undirected edge matches either orientation); EdgeDelete[g, patt] removes every edge matching patt. Vertices, order and the remaining weights are kept.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= EdgeList[EdgeDelete[CycleGraph[4], 1<->2]]
Out[1]= {2 <-> 3, 3 <-> 4, 4 <-> 1}

In[2]:= EdgeList[EdgeDelete[CycleGraph[4], {2<->1, 3<->4}]]
Out[2]= {2 <-> 3, 4 <-> 1}

In[3]:= EdgeList[EdgeDelete[CycleGraph[4], _[1, _]]]
Out[3]= {2 <-> 3, 3 <-> 4, 4 <-> 1}

In[4]:= EdgeDelete[CycleGraph[4], 1->2]
Out[4]= EdgeDelete[Graph[<4 vertices, 4 edges>], 1 -> 2]

In[5]:= EdgeDelete[CycleGraph[4], 1<->3]
Out[5]= EdgeDelete[Graph[<4 vertices, 4 edges>], TwoWayRule[1, 3]]
```

## Implementation notes

- `Protected`. A non-graph first argument is left unevaluated.
- Each listed edge must exist, otherwise the call is left unevaluated. An
  undirected edge matches either orientation; `1 -> 2` is not an edge of an
  undirected graph.
- Orders and the remaining weights are kept.
- `O(V + E)` plus one hash per argument item; result seeded into the graph memo.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
