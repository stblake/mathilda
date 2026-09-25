# PlanarGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PlanarGraphQ[g] gives True if g can be drawn in the plane without edge crossings. Linear-time left-right planarity test; gives False for non-graphs.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= PlanarGraphQ[CompleteGraph[4]]
Out[1]= True

In[2]:= PlanarGraphQ[CompleteGraph[5]]
Out[2]= False

In[3]:= PlanarGraphQ[CompleteGraph[{3,3}]]
Out[3]= False

In[4]:= PlanarGraphQ[PetersenGraph[]]
Out[4]= False

In[5]:= PlanarGraphQ[Graph[{1->2,2->3,3->1}]]
Out[5]= True

In[6]:= PlanarGraphQ[x]
Out[6]= False
```

## Implementation notes

- `Protected`; gives `False` for a non-graph. Edge directions are ignored.
- Algorithm (`galg_planar.c`): the linear-time left-right planarity test of de
  Fraysseix and Rosenstiehl in Brandes' formulation (the same algorithm as
  networkx's `check_planarity`), without the embedding phase. Two iterative DFS
  passes per connected component (orientation with lowpoints and nesting
  depths, then conflict-pair testing), so path-like graphs with 10^6 vertices
  need no call stack. The Euler bound rejects `m > 3n - 6` up front.
- O(n + m) time and memory.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
