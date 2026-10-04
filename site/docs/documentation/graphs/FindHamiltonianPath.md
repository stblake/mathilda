# FindHamiltonianPath

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindHamiltonianPath[g] gives a Hamiltonian path of g as a vertex list, or {} if there is none; FindHamiltonianPath[g, s, t] one from s to t. Exact; unevaluated if the search budget is exhausted.`**

## Examples (9)

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

### Applications (3)

A cycle gives a path through every vertex

```mathematica
In[7]:= FindHamiltonianPath[CycleGraph[5]]
Out[7]= {5, 1, 2, 3, 4}
```

The path itself

```mathematica
In[8]:= FindHamiltonianPath[PathGraph[{1, 2, 3}]]
Out[8]= {1, 2, 3}
```

A Hamiltonian path from 1 to 3

```mathematica
In[9]:= FindHamiltonianPath[CompleteGraph[4], 1, 3]
Out[9]= {1, 2, 4, 3}
```

## Implementation notes

**Algorithm.** `builtin_find_hamiltonian_path` finds a path visiting every vertex
exactly once, as a vertex list, or `{}` if none exists; `FindHamiltonianPath[g,
s, t]` one running from `s` to `t`. It reduces to a Hamiltonian *cycle* search by
adding a virtual vertex `z`: joined to every vertex for the unconstrained form,
or only to `s` and `t` (as `z -> s` and `t -> z`, directed-aware) for the fixed-
endpoint form, so a cycle through `z` is a path between the required ends. The
cycle search itself branches over **edge decisions** with constraint propagation:
each vertex needs exactly two chosen incident edges (one in-arc and one out-arc
when directed), a chosen edge closing a fragment shorter than `n` is refused, and
at every node the remaining non-excluded graph must stay biconnected (strongly
connected when directed, Tarjan, linear). Every step is trailed and undone, so
`{}` is a proof.

**Data structures.** A `Ham` instance with CSR incidence lists (out/in arcs when
directed), per-vertex chosen/available counts, a fragment end-link array
`oth[]`, a propagation queue, Tarjan scratch, and an explicit decision trail
(`HamTr`). Items come from the memo's integer endpoints (`graph_edge_indices`),
deduplicated.

**Complexity / limits.** Easy on small sparse instances (a random cubic
400-vertex graph ~1 ms); the per-node biconnectivity check is linear, so dense
meshes grow roughly quadratically. A budget of `HAM_MAX_NODES` search nodes,
polled against the `TimeConstrained` deadline, leaves the call unevaluated when
exhausted — never a wrong "no".

- `Protected`; unevaluated on a non-graph.
- The one-vertex graph gives `{}`, as in Mathematica.
- Paths reduce to cycles through an added vertex and use the
  `FindHamiltonianCycle` engine (exact; polls `TimeConstrained`).

**Attributes:** `Protected`.

## References

**See also:** [FindHamiltonianCycle](../../graphs/FindHamiltonianCycle/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- Source: [`src/graph/galg_hamilton.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_hamilton.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

The result is a list of vertices in visiting order, or `{}` when the graph has
no Hamiltonian path. `FindHamiltonianPath[g, s, t]` additionally fixes the two
endpoints.

The search is exact and complete, so `{}` is a proof of non-existence — but it
is backtracking, so a hard instance can exhaust the search budget and leave the
call unevaluated rather than return a wrong answer; bound it with
`TimeConstrained`.
