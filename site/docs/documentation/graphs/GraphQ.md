# GraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphQ[g] gives True if g is a valid graph, and False otherwise.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= GraphQ[Graph[{1,2}, {1->2}]]
Out[1]= True

In[2]:= GraphQ[Graph[{1}, {1->1}]]
Out[2]= False

In[3]:= GraphQ[5]
Out[3]= False

In[4]:= GraphQ[Graph[{}, {}]]
Out[4]= True

In[5]:= GraphQ[CycleGraph[4]]
Out[5]= True
```

## Algorithm

graphq.c - GraphQ[g]: is g a valid graph?

A thin wrapper over graph_is_valid (graph_util.c): returns the symbol True when the (already-evaluated) argument is a canonical, valid graph, and False otherwise. Non-unary calls are left unevaluated (NULL).

Memory (SPEC section 4): returns a freshly-allocated symbol; the evaluator frees `res`.

## Implementation notes

- `Protected`. Never left unevaluated: any non-graph gives `False`.
- A graph is valid when it is the canonical `Graph[List, List]` with every edge
  a 2-argument `DirectedEdge`/`UndirectedEdge`, no self-loops, no parallel
  edges, and every endpoint present in the vertex list (the same conditions the
  `Graph` constructor enforces).
- The null graph `Graph[{}, {}]` is valid.

**Attributes:** `Protected`.

## References

**See also:** [Graph](../../graphs/Graph/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)
