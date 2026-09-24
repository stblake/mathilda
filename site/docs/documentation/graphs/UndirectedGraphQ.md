# UndirectedGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`UndirectedGraphQ[g] gives True if every edge of g is undirected (including when g has no edges), and False otherwise.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= UndirectedGraphQ[Graph[{1,2},{}]]
Out[1]= True

In[2]:= UndirectedGraphQ[CycleGraph[3]]
Out[2]= True

In[3]:= UndirectedGraphQ[Graph[{1,2,3},{1->2,2<->3}]]
Out[3]= False

In[4]:= UndirectedGraphQ[x]
Out[4]= False
```

## Implementation notes

- `Protected`. Edgeless graphs are undirected; a mixed graph is neither directed
  nor undirected (see `DirectedGraphQ`).
- Shared by all the structural predicates (`UndirectedGraphQ`, `EmptyGraphQ`,
  `CompleteGraphQ`, `BipartiteGraphQ`, `VertexQ`, `EdgeQ`, `AcyclicGraphQ`,
  `TreeGraphQ`): each gives `False` — never unevaluated — for an argument that
  is not a valid graph, and each runs in linear time on first use and `O(1)`
  when repeated on the same graph (see the *Performance model* section of this
  file's preamble).

**Attributes:** `Protected`.

## References

**See also:** [DirectedGraphQ](../../graphs/DirectedGraphQ/), [EmptyGraphQ](../../graphs/EmptyGraphQ/), [CompleteGraphQ](../../graphs/CompleteGraphQ/), [BipartiteGraphQ](../../graphs/BipartiteGraphQ/), [VertexQ](../../graphs/VertexQ/), [EdgeQ](../../graphs/EdgeQ/), [AcyclicGraphQ](../../graphs/AcyclicGraphQ/), [TreeGraphQ](../../graphs/TreeGraphQ/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
