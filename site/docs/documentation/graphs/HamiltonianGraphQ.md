# HamiltonianGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HamiltonianGraphQ[g] gives True if g has a Hamiltonian cycle. Exact; unevaluated if the search budget is exhausted.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= HamiltonianGraphQ[CompleteGraph[5]]
Out[1]= True

In[2]:= HamiltonianGraphQ[PetersenGraph[]]
Out[2]= False

In[3]:= HamiltonianGraphQ[Graph[{1},{}]]
Out[3]= True

In[4]:= HamiltonianGraphQ[Graph[{},{}]]
Out[4]= False

In[5]:= HamiltonianGraphQ[x]
Out[5]= False
```

## Implementation notes

- `Protected`; gives `False` for a non-graph.
- `True` for the one-vertex graph, `False` for the graph with no vertices.
- Uses the `FindHamiltonianCycle` constraint-propagation search (exact;
  polls `TimeConstrained`).

**Attributes:** `Protected`.

## References

**See also:** [FindHamiltonianCycle](../../graphs/FindHamiltonianCycle/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
