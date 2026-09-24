# EdgeIndex

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeIndex[g, e] gives the position of the edge e in EdgeList[g] (an undirected edge matches either orientation; u->v and u<->v sugar are accepted); EdgeIndex[g, {e1, ...}] gives a list of positions.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= VertexIndex[Graph[{a<->b, b<->c}], b]
Out[1]= 2

In[2]:= VertexIndex[Graph[{a<->b, b<->c}], {c, a}]
Out[2]= {3, 1}

In[3]:= EdgeIndex[CycleGraph[4], 3<->2]
Out[3]= 2

In[4]:= EdgeIndex[CycleGraph[4], 3->4]
Out[4]= 3

In[5]:= EdgeIndex[CycleGraph[4], {1<->2, 4<->1}]
Out[5]= {1, 4}

In[6]:= EdgeIndex[Graph[{1->2}], 2->1]
Out[6]= EdgeIndex[Graph[<2 vertices, 1 edge>], 2 -> 1]
```

## Implementation notes

- `Protected`. A non-graph first argument, or a vertex/edge not in the graph, is
  left unevaluated.
- `EdgeIndex` matches `u -> v` against an undirected edge of an undirected
  graph, as Mathematica does; an undirected edge matches either orientation.

**Attributes:** `Protected`.

## References

**See also:** [VertexIndex](../../graphs/VertexIndex/)

- Source: [`src/graph/gops_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gops_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
