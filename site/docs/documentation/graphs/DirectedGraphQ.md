# DirectedGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DirectedGraphQ[g] gives True if g has at least one edge and all edges of g are directed. An edgeless graph counts as undirected.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= DirectedGraphQ[Graph[{1->2,2->3}]]
Out[1]= True

In[2]:= DirectedGraphQ[CycleGraph[3]]
Out[2]= False

In[3]:= DirectedGraphQ[Graph[{1,2},{}]]
Out[3]= False

In[4]:= DirectedGraphQ[Graph[{1,2,3},{1->2,2<->3}]]
Out[4]= False

In[5]:= DirectedGraphQ[5]
Out[5]= False
```

## Algorithm

directedq.c - DirectedGraphQ[g]: True iff g is a valid graph with at least one edge, all of them DirectedEdge. False otherwise -- including for an edgeless graph, which (as in the Wolfram Language) counts as undirected, so DirectedGraphQ and UndirectedGraphQ are never both True. Memory (SPEC section 4): returns a fresh symbol; the evaluator frees res.

## Implementation notes

- `Protected`. `False` (never unevaluated) for a non-graph.
- An edgeless graph counts as undirected (as in the Wolfram Language), so
  `DirectedGraphQ` and `UndirectedGraphQ` are never both `True`. A mixed graph
  is neither.

**Attributes:** `Protected`.

## References

**See also:** [UndirectedGraphQ](../../graphs/UndirectedGraphQ/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
