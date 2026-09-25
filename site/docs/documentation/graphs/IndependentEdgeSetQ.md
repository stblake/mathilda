# IndependentEdgeSetQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IndependentEdgeSetQ[g, elist] gives True if elist is a set of edges of g no two of which share a vertex.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= IndependentEdgeSetQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[4,3]}]
Out[1]= True

In[2]:= IndependentEdgeSetQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[2,3]}]
Out[2]= False

In[3]:= IndependentEdgeSetQ[Graph[{1->2,3->4}], {DirectedEdge[2,1]}]
Out[3]= False

In[4]:= EdgeCoverQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[3,4]}]
Out[4]= True

In[5]:= EdgeCoverQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[2,3]}]
Out[5]= False

In[6]:= EdgeCoverQ[CycleGraph[4], {UndirectedEdge[1,3]}]
Out[6]= False
```

## Implementation notes

- `Protected` membership predicates; give `False` when `g` is not a graph.
- An element that is not an edge of `g` gives `False`. An `UndirectedEdge`
  matches either orientation; a `DirectedEdge` matches only as given.

**Attributes:** `Protected`.

## References

**See also:** [EdgeCoverQ](../../graphs/EdgeCoverQ/)

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
