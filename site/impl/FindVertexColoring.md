---
references:
  - "D. Brélaz, *New methods to color the vertices of a graph*, Communications of the ACM **22** (1979) 251-256."
source: src/graph/vertexcoloring.c
---
**Algorithm.** `builtin_find_vertex_coloring` returns a colour assignment whose
number of distinct colours equals the **chromatic number** — Wolfram's documented
meaning, and the reason the search is exact rather than greedy: a merely-valid
colouring with too many colours would be a plausible list of integers that
silently contradicts that contract. It brackets the chromatic number between a
cheap upper bound `ub` from a **DSATUR** pass (`fvc_dsatur_bound`, which also
exhibits a real colouring) and a lower bound `lb` from a multi-start greedy clique
(`fvc_clique_bound`). When `lb == ub` the DSATUR colouring is proven optimal and
returned with zero search — this is what makes `CompleteGraph[128]` immediate
rather than a hang. Otherwise it runs a **DSATUR branch-and-bound** (`fvc_bb`) that
picks the next vertex dynamically as the uncoloured one of maximum saturation and
prunes on three grounds: bound (a partial colouring already at `best_k` colours
cannot win), symmetry breaking (a vertex tries only colours `1..used+1`), and
optimality (stop once `best_k == lb`). Adjacency is the **undirected**
neighbourhood — an edge constrains its endpoints whichever way it points. The empty
graph colours to `{}`. Exactness is preserved by **refusing** (returning `NULL`,
leaving the call unevaluated) whenever minimality cannot be proven: above
`FVC_MAX_VERTICES = 128`, when the node budget `FVC_MAX_STEPS = 8,000,000` is
spent, or on allocation failure — never by returning the incumbent. A second
argument also declines (the `FindVertexColoring[g, {c1, ...}]` forms are a later
layer). This is Wolfram's own `"BacktrackingDS"` method, so shipping only it is a
documented subset.

**Data structures.** `graph_build_adj` yields a `GraphAdj` with successor `out[]`
and predecessor `in[]` lists (an undirected edge appears in both), walked in place
so there is nothing beyond the `GraphAdj` to free. The working colouring `col[]`
is 1-based (0 = uncoloured); the branch-and-bound threads an `FvcBB` struct
(working colouring, best complete colouring, incumbent `best_k`, lower bound `lb`)
by pointer, with `seen[]`/`forbid[]` as stack arrays sized to the vertex cap so the
hot path never allocates. The result is a `List` of positive-integer colour labels,
one per vertex in `VertexList` order.

**Complexity / limits.** Exact graph colouring is NP-hard; the vertex cap bounds
size and the **node count** bounds cost *deterministically* — a wall-clock cutoff
would make the answer machine-dependent (a fast host proves minimality, a slow one
refuses the same graph), whereas a fixed node budget gives every machine the same
answer. `fvc_bb` polls `tc_check_deadline()` every 4096 nodes so an interactive
`TimeConstrained[...]` is honoured even where `SIGPROF` is unreliable; the node
budget is only the backstop for an unattended run.
