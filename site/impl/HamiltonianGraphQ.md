---
source: src/graph/galg_hamilton.c
---
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
