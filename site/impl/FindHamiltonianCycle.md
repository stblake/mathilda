---
references:
  - "R. E. Tarjan, *Depth-first search and linear graph algorithms*, SIAM J. Comput. **1** (1972) 146-160."
source: src/graph/galg_hamilton.c
---
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
