# FindVertexCover

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindVertexCover[g] gives a minimum vertex cover of g (a smallest vertex set touching every edge), exact by branch and reduce. Stays unevaluated if optimality cannot be proven within the search budget.`**

## Examples (10)

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

### Applications (6)

The centre alone covers every spoke

```mathematica
In[5]:= FindVertexCover[StarGraph[5]]
Out[5]= {1}
```

Works on arbitrary vertex names

```mathematica
In[6]:= FindVertexCover[PathGraph[{a, b, c, d, e}]]
Out[6]= {b, d}
```

All but one vertex

```mathematica
In[7]:= FindVertexCover[CompleteGraph[4]]
Out[7]= {1, 2, 3}
```

Odd cycle needs three

```mathematica
In[8]:= FindVertexCover[CycleGraph[5]]
Out[8]= {1, 2, 4}
```

The result passes the predicate

```mathematica
In[9]:= VertexCoverQ[CycleGraph[5], FindVertexCover[CycleGraph[5]]]
Out[9]= True
```

Minimum cover has size ten minus alpha

```mathematica
In[10]:= Length[FindVertexCover[PetersenGraph[]]]
Out[10]= 6
```

## Implementation notes

**Algorithm.** `builtin_find_vertex_cover` takes exactly one argument and returns a *minimum* vertex cover. A set is a vertex cover exactly when its complement is independent, so it runs the same exact solver as `FindIndependentVertexSet` (`galg_max_independent_set`) and returns the vertices *not* in the maximum independent set it finds. Edge direction is ignored. The solver is per component: König's theorem for large bipartite components, the complement-clique bound for dense ones, and branch and reduce (degree-0/1/2, folding, domination, clique-cover bound, mirror branching) for sparse ones.

**Data structures.** The `Graph[List, List]` expression is flattened to a `GalgUG` CSR adjacency; a byte mask marks the independent set and the cover is read off as the unmarked vertices, in vertex order, and rebuilt as a list of the graph's own vertex expressions.

**Complexity / limits.** Minimum vertex cover is NP-hard, so the worst case is exponential. The solver gives up after a fixed node/work/memory/depth budget or when the `TimeConstrained` deadline passes, and the head then stays unevaluated; it never returns a non-minimum cover. A graph that is not valid, or any form with more than one argument, is left unevaluated. Ties between optimal covers are broken deterministically by the vertex order.

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

- F. V. Fomin, F. Grandoni and D. Kratsch, *A measure and conquer approach for the analysis of exact algorithms*, J. ACM **56**(5) (2009) Art. 25.
- Source: [`src/graph/galg_mis.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_mis.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

The cover is a minimum one: it is the complement of a maximum independent set, found by the same exact solver as `FindIndependentVertexSet`. When several minimum covers exist the one returned is fixed by the vertex order, so repeated calls agree, but it need not be the one another system would pick.

Edge direction is ignored. Minimum vertex cover is NP-hard; if the search budget or a `TimeConstrained` limit is exhausted the call stays unevaluated rather than returning a non-minimum cover.
