# BipartiteGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BipartiteGraphQ[g] gives True if the vertices of g split into two sets with every edge running between them (edge direction is ignored).`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= BipartiteGraphQ[CycleGraph[5]]
Out[1]= False

In[2]:= BipartiteGraphQ[CycleGraph[6]]
Out[2]= True

In[3]:= BipartiteGraphQ[Graph[{1,2,3},{}]]
Out[3]= True

In[4]:= BipartiteGraphQ[x]
Out[4]= False
```

## Implementation notes

- `Protected`. Edge direction is ignored; edgeless graphs are bipartite.
  `False` for a non-graph (see `UndirectedGraphQ`).

**Attributes:** `Protected`.

## References

**See also:** [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)
