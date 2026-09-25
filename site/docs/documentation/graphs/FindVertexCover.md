# FindVertexCover

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindVertexCover[g] gives a minimum vertex cover of g (a smallest vertex set touching every edge), exact by branch and reduce. Stays unevaluated if optimality cannot be proven within the search budget.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= FindVertexCover[PetersenGraph[]]
Out[1]= {1, 4, 5, 7, 8, 10}

In[2]:= FindVertexCover[StarGraph[5]]
Out[2]= {1}

In[3]:= FindVertexCover[CycleGraph[5]]
Out[3]= {1, 2, 4}

In[4]:= FindVertexCover[x]
Out[4]= FindVertexCover[x]
```

## Implementation notes

- `Protected`; unevaluated on a non-graph.
- Computed as the complement of a maximum independent set, so it shares the
  exact independent-set engine and its budgets (see `FindIndependentVertexSet`).
- **Exactness policy** (all NP-hard heads: `FindVertexCover`,
  `FindIndependentVertexSet`, `FindClique`, `FindKClique`, the Hamiltonian heads,
  and the isomorphism search): the result is a *proven* optimum / a complete
  answer, or the head stays unevaluated when its deterministic node budget is
  exhausted. All of them poll the `TimeConstrained` deadline, so e.g.
  `TimeConstrained[FindClique[g], 1]` returns `$Aborted` rather than hanging.
  They never return a merely-good answer.
- Weighted graphs: optimizes cardinality, as the Wolfram documentation states.

**Attributes:** `Protected`.

## References

**See also:** [VertexList](../../graphs/VertexList/), [FindIndependentVertexSet](../../graphs/FindIndependentVertexSet/), [FindClique](../../graphs/FindClique/), [FindKClique](../../graphs/FindKClique/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
