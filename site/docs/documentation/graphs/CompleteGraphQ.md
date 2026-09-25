# CompleteGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CompleteGraphQ[g] gives True if every pair of distinct vertices of g is joined by an edge in both directions (an undirected edge, or directed edges both ways). CompleteGraphQ[g, vlist] tests the subgraph induced by vlist; False if some element of vlist is not a vertex of g.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= CompleteGraphQ[Graph[{1->2,2->1}]]
Out[1]= True

In[2]:= CompleteGraphQ[Graph[{1->2}]]
Out[2]= False

In[3]:= CompleteGraphQ[CompleteGraph[5]]
Out[3]= True

In[4]:= CompleteGraphQ[CycleGraph[4], {1,2}]
Out[4]= True

In[5]:= CompleteGraphQ[CycleGraph[4], {1,3}]
Out[5]= False

In[6]:= CompleteGraphQ[CycleGraph[4], {1,9}]
Out[6]= False
```

## Implementation notes

- `Protected`. A pair `(u, v)` is covered by an undirected edge or a directed
  `u -> v`, so a complete directed graph needs both directions.
- Graphs with 0 or 1 vertices are complete.
- `CompleteGraphQ[g, vlist]` gives `False` if some element is not a vertex of
  `g`; repeats are ignored; `{}` is complete.
- `False` for a non-graph (see `UndirectedGraphQ`).

**Attributes:** `Protected`.

## References

**See also:** [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
