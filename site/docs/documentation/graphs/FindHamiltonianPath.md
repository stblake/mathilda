# FindHamiltonianPath

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindHamiltonianPath[g] gives a Hamiltonian path of g as a vertex list, or {} if there is none; FindHamiltonianPath[g, s, t] one from s to t. Exact; unevaluated if the search budget is exhausted.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= FindHamiltonianPath[PetersenGraph[]]
Out[1]= {6, 7, 2, 4, 1, 3, 5, 10, 9, 8}

In[2]:= FindHamiltonianPath[CycleGraph[5], 1, 2]
Out[2]= {1, 5, 4, 3, 2}

In[3]:= FindHamiltonianPath[CycleGraph[5], 1, 3]
Out[3]= {}

In[4]:= FindHamiltonianPath[Graph[{1->2,2->3}]]
Out[4]= {1, 2, 3}

In[5]:= FindHamiltonianPath[StarGraph[4]]
Out[5]= {}

In[6]:= FindHamiltonianPath[Graph[{1},{}]]
Out[6]= {}
```

## Implementation notes

- `Protected`; unevaluated on a non-graph.
- The one-vertex graph gives `{}`, as in Mathematica.
- Paths reduce to cycles through an added vertex and use the
  `FindHamiltonianCycle` engine (exact; polls `TimeConstrained`).

**Attributes:** `Protected`.

## References

**See also:** [FindHamiltonianCycle](../../graphs/FindHamiltonianCycle/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
