# HamiltonianGraphQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HamiltonianGraphQ[g] gives True if g has a Hamiltonian cycle. Exact; unevaluated if the search budget is exhausted.`**

## Examples (9)

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

### Applications (4)

A cycle is its own Hamiltonian cycle

```mathematica
In[6]:= HamiltonianGraphQ[CycleGraph[5]]
Out[6]= True
```

Every complete graph on >= 3 vertices is Hamiltonian

```mathematica
In[7]:= HamiltonianGraphQ[CompleteGraph[4]]
Out[7]= True
```

A path has no cycle at all

```mathematica
In[8]:= HamiltonianGraphQ[PathGraph[{1, 2, 3}]]
Out[8]= False
```

The classic non-Hamiltonian example

```mathematica
In[9]:= HamiltonianGraphQ[PetersenGraph[]]
Out[9]= False
```

## Implementation notes

**Algorithm.** `builtin_hamiltonian_graph_q` is `True` iff `g` has a Hamiltonian
cycle. It asks the shared cycle search (`ham_cycles`) for a single cycle and
reports whether one was found. The search branches over **edge decisions** with
constraint propagation — each vertex needs exactly two chosen incident edges, an
edge closing a fragment into a cycle shorter than `n` is refused, and the
remaining non-excluded graph must stay biconnected (strongly connected when
directed) at every node (Tarjan, linear) — with every step trailed and undone, so
a `False` answer is a completeness proof, not a heuristic. Small cases are
special-cased: the one-vertex graph is Hamiltonian (empty cycle), `K2` is not,
and a directed 2-cycle is.

**Data structures.** The `Ham` instance: CSR incidence (out/in arcs when
directed), chosen/available per-vertex counters, a fragment end-link array, a
propagation queue, Tarjan scratch, and a decision trail. Edge items come from the
memo's integer endpoints, deduplicated.

**Complexity / limits.** Easy on small sparse graphs; the linear per-node
connectivity check makes dense meshes grow roughly quadratically. A budget of
`HAM_MAX_NODES` nodes, polled against the `TimeConstrained` deadline, leaves the
call **unevaluated** when exhausted rather than return a wrong `False`. A
non-graph argument gives `False`.

- `Protected`; gives `False` for a non-graph.
- `True` for the one-vertex graph, `False` for the graph with no vertices.
- Uses the `FindHamiltonianCycle` constraint-propagation search (exact;
  polls `TimeConstrained`).

**Attributes:** `Protected`.

## References

**See also:** [FindHamiltonianCycle](../../graphs/FindHamiltonianCycle/), [TimeConstrained](../../time-and-date/TimeConstrained/)

- Source: [`src/graph/galg_hamilton.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_hamilton.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

`HamiltonianGraphQ[g]` tests for a Hamiltonian *cycle* (every vertex visited
once, returning to the start). The Petersen graph is the standard example of a
graph that is 3-regular and connected yet non-Hamiltonian.

The decision is exact, so `False` is a proof — but it is a backtracking search,
so a hard instance can exhaust the budget and leave the call unevaluated rather
than return a wrong answer; bound it with `TimeConstrained`.
