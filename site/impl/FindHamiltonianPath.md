---
source: src/graph/galg_hamilton.c
---
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
