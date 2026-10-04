# FindHamiltonianCycle

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindHamiltonianCycle[g] gives {c} with c a Hamiltonian cycle of g as a list of edges, or {} if there is none. FindHamiltonianCycle[g, n] / [g, All] gives up to n / all Hamiltonian cycles. Exact backtracking with pruning; unevaluated if the search budget is exhausted.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= FindHamiltonianCycle[CycleGraph[4]]
Out[1]= {{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}}

In[2]:= FindHamiltonianCycle[CompleteGraph[4], All]
Out[2]= {{1 <-> 2, 2 <-> 4, 4 <-> 3, 3 <-> 1}, {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}, {1 <-> 3, 3 <-> 2, 2 <-> 4, 4 <-> 1}}

In[3]:= FindHamiltonianCycle[CompleteGraph[5], 2]
Out[3]= {{1 <-> 2, 2 <-> 4, 4 <-> 5, 5 <-> 3, 3 <-> 1}, {1 <-> 2, 2 <-> 5, 5 <-> 4, 4 <-> 3, 3 <-> 1}}

In[4]:= FindHamiltonianCycle[PetersenGraph[]]
Out[4]= {}

In[5]:= FindHamiltonianCycle[Graph[{1->2,2->1}]]
Out[5]= {{1 -> 2, 2 -> 1}}

In[6]:= FindHamiltonianCycle[Graph[{1},{}]]
Out[6]= {{}}
```

### Scope (2)

```mathematica
In[7]:= FindHamiltonianCycle[CompleteGraph[2]]
Out[7]= {}

In[8]:= TimeConstrained[FindHamiltonianCycle[GridGraph[{150, 150}]], 0.05]
Out[8]= $Aborted
```

### Applications (2)

The cycle itself, as a list of edges

```mathematica
In[9]:= FindHamiltonianCycle[CycleGraph[5]]
Out[9]= {{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 5, 5 <-> 1}}
```

A path has no Hamiltonian cycle

```mathematica
In[10]:= FindHamiltonianCycle[PathGraph[4]]
Out[10]= {}
```

## Implementation notes

**Algorithm.** `builtin_find_hamiltonian_cycle` is an exact backtracking search over *edge*
decisions with constraint propagation. Each edge is undecided, chosen or excluded on an undo
trail. Propagation enforces that every vertex has exactly two chosen incident edges (for a
directed graph, one in-arc and one out-arc): once two are chosen the rest are excluded, and when
only two remain undecided both are forced; fewer than two is a contradiction. Chosen edges form
vertex-disjoint path fragments, and an edge that would close a fragment into a sub-tour of fewer
than `n` vertices is refused. The strong pruning is a connectivity test at every node — the
non-excluded graph must stay biconnected (undirected, iterative Tarjan) or strongly connected
(directed, double BFS); branching extends the fragment end with the fewest remaining options.

**Data structures.** A CSR incidence (`off`/`inc`, plus `ioff`/`iinc` for directed graphs);
an edge-state array `st[]`; per-vertex chosen/available counts; a fragment-end map `oth[]`; an
undo trail; a propagation worklist; and Tarjan scratch arrays (`disc`/`low`/`stk`).

**Complexity / limits.** NP-complete; the per-node biconnectivity check is linear. Capped at
`HAM_MAX_NODES = 2·10^7` search nodes, polled against `TimeConstrained`; exhaustion returns
unevaluated (never a wrong "no cycle"). Returns `{c}` with `c` the cycle as a list of edges
(`FindHamiltonianCycle[g, n]` / `[g, All]` returns up to `n` / all cycles), or `{}` when none
exists. A non-graph argument returns unevaluated.

- `Protected`; unevaluated on a non-graph.
- The cycle starts at `VertexList[g][[1]]`, each edge written in traversal
  order.
- The one-vertex graph gives `{{}}`; `K2` has none; a directed 2-cycle is one.
- Exact (see the exactness policy under `FindVertexCover`): complete answer or
  unevaluated on budget exhaustion; polls `TimeConstrained`.
- Algorithm: search over edge decisions with constraint propagation: each vertex
  needs exactly two chosen edges (directed: one in, one out), chosen edges form
  path fragments whose end-to-end links forbid short cycles, and the remaining
  graph must stay biconnected (directed: strongly connected) at every node.
- Performance: a random 400-vertex cubic graph takes about 1 ms. Large meshes
  are not a fast case: the linear per-node biconnectivity check makes
  `GridGraph[{n, n}]` roughly quadratic in the vertex count (100 x 100 about
  1.5 s, 150 x 150 about 7 s), so wrap far larger instances in
  `TimeConstrained`.

**Attributes:** `Protected`.

## References

**See also:** [FindVertexCover](../../graphs/FindVertexCover/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- R. E. Tarjan, *Depth-first search and linear graph algorithms*, SIAM J. Comput. **1** (1972) 146-160.
- Source: [`src/graph/galg_hamilton.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_hamilton.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

The result is `{c}`, where `c` is a Hamiltonian cycle — a closed tour visiting every vertex
once — given as a list of edges in traversal order; it is `{}` when no such cycle exists.
`FindHamiltonianCycle[g, n]` returns up to `n` cycles and `[g, All]` returns all of them.

The search is exact backtracking with biconnectivity/strong-connectivity pruning, so both a
cycle and the answer `{}` are proofs. If the node budget is exhausted the call stays
unevaluated.
