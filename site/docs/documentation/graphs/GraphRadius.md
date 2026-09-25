# GraphRadius

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphRadius[g] gives the minimum vertex eccentricity of g. Infinity unless g is (strongly) connected.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= GraphDiameter[PathGraph[5]]
Out[1]= 4

In[2]:= GraphRadius[PathGraph[5]]
Out[2]= 2

In[3]:= GraphDiameter[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= Infinity

In[4]:= GraphDiameter[Graph[{},{}]]
Out[4]= 0
```

### Options (2)

```mathematica
In[5]:= GraphRadius[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}]]
Out[5]= 3.0

In[6]:= GraphDiameter[Graph[{1,2,3},{1<->2,2<->3,1<->3}, EdgeWeight->{1,1,5}]]
Out[6]= 2.0
```

## Implementation notes

- *(w)* weight-aware (machine reals when weighted).
- Unweighted and not strongly connected (not connected, if undirected):
  `Infinity`.
- Weighted: taken over the weighted eccentricities, so a weighted digraph that
  is not strongly connected can have a finite radius
  (*reverse-engineered*).
- Graphs with no vertices give diameter/radius `0`.
- Reduced from the cached MS-BFS per-source summary, so calling several of
  `GraphDiameter`, `GraphRadius`, `GraphCenter`, `GraphPeriphery`,
  `ClosenessCentrality`, `EccentricityCentrality` on one graph pays for one
  all-pairs pass.

**Attributes:** `Protected`.

## References

**See also:** [GraphDiameter](../../graphs/GraphDiameter/), [GraphCenter](../../graphs/GraphCenter/), [GraphPeriphery](../../graphs/GraphPeriphery/), [ClosenessCentrality](../../graphs/ClosenessCentrality/), [EccentricityCentrality](../../graphs/EccentricityCentrality/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
