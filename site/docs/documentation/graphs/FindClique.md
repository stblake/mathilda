# FindClique

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindClique[g] gives {c} with c a maximum clique of g. FindClique[g, k] finds a maximal clique of at most k vertices, [g, {k}] of exactly k, [g, {kmin, kmax}] within the range; a third argument n (or All) gives up to n such cliques. For directed graphs a clique needs edges both ways. Exact branch and bound with colouring bounds.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= FindClique[CompleteGraph[5]]
Out[1]= {{1, 2, 3, 4, 5}}

In[2]:= FindClique[CycleGraph[5], Infinity, All]
Out[2]= {{4, 5}, {3, 4}, {2, 3}, {1, 5}, {1, 2}}

In[3]:= FindClique[WheelGraph[6], {3}, 2]
Out[3]= {{1, 2, 6}, {1, 2, 3}}

In[4]:= FindClique[Graph[{1,2,3,4},{UndirectedEdge[1,2],UndirectedEdge[2,3],UndirectedEdge[1,3],UndirectedEdge[3,4]}], {2,3}, All]
Out[4]= {{1, 2, 3}, {3, 4}}

In[5]:= FindClique[CompleteGraph[4], 2]
Out[5]= {}

In[6]:= FindClique[Graph[{1->2,2->1,2->3,3->2,1->3}]]
Out[6]= {{1, 2}}
```

### Scope (1)

```mathematica
In[7]:= TimeConstrained[FindClique[RandomGraph[{400, 40000}]], 0.001]
Out[7]= $Aborted
```

## Options & behaviour

A search that exceeds its `TimeConstrained` limit aborts instead of hanging:

## Implementation notes

- `Protected`; unevaluated on a non-graph.
- A directed graph needs edges both ways between clique members.
- Maximal cliques are listed by size, largest first; within one size
  Mathematica's order is reproduced
  (`FindClique[CycleGraph[5], Infinity, All]` gives
  `{{4,5},{3,4},{2,3},{1,5},{1,2}}`).
- `FindClique[g, 2]` is `{}` when every maximal clique is larger.
- Exact (see the exactness policy under `FindVertexCover`): a proven maximum or
  unevaluated when the deterministic node budget runs out; the
  `TimeConstrained` deadline is polled, so a timed-out search gives `$Aborted`.
- Maximum clique: bitset branch and bound with greedy-colouring bounds (Tomita's
  MCQ/MCS family in San Segundo's BBMC form), vertices numbered by reverse
  degeneracy order; graphs above 3000 vertices are decomposed by degeneracy
  (each vertex's later neighbourhood is a small local problem, skipped when its
  core number cannot beat the incumbent).
- Enumeration: pivoted Bron-Kerbosch on bitsets with size pruning.

**Attributes:** `Protected`.

## References

**See also:** [FindVertexCover](../../graphs/FindVertexCover/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
