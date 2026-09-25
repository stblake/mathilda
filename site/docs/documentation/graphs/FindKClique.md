# FindKClique

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindKClique[g, k] gives {c} with c a largest k-clique of g: a vertex set whose members are pairwise within distance k.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= FindKClique[CycleGraph[8], 2]
Out[1]= {{1, 2, 3}}

In[2]:= FindKClique[PathGraph[{1,2,3,4,5}], 1]
Out[2]= {{1, 2}}

In[3]:= FindKClique[x, 2]
Out[3]= FindKClique[x, 2]
```

## Implementation notes

- `Protected`; unevaluated on a non-graph.
- Computed as a maximum clique of the `k`-th power graph, using the same exact
  BBMC maximum-clique search as `FindClique` (proven optimum or unevaluated on
  budget exhaustion; polls `TimeConstrained`).

**Attributes:** `Protected`.

## References

**See also:** [FindClique](../../graphs/FindClique/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
