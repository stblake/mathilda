---
references:
  - "A. B. Kahn, *Topological sorting of large networks*, Comm. ACM **5** (1962) 558-562."
source: src/graph/acyclic.c
---
**Algorithm.** `builtin_acyclic_graph_q` decides whether the graph has a cycle, treating
directed edges as forward-only and undirected edges as traversable either way, by an exact
three-step contraction. (1) A **union-find** (path halving, union by size) joins the endpoints
of every undirected edge; a join of two already-connected vertices closes an undirected cycle,
so the graph is cyclic. (2) Each directed edge is mapped onto the contracted components; a
directed edge inside one component also closes a cycle. (3) The contracted digraph is run
through **Kahn's topological sort** — the graph is acyclic iff every vertex is eventually
removed. So an undirected graph is acyclic exactly when it is a forest, and a directed graph
exactly when it is a DAG; an anti-parallel pair `u -> v, v -> u` is a 2-cycle. The 0/1 answer is
cached on the graph node (`GRAPH_PROP_ACYCLIC`).

**Data structures.** The union-find `parent[]`/`size[]`; pre-resolved endpoint index arrays
(`eu`/`ev`/`edir`) from the validated-graph memo; staged directed tails/heads and in/out degree
arrays; a CSR `start[]`/`succ[]` for the contracted digraph and a Kahn work queue.

**Complexity / limits.** `O(V + E·alpha(V))`, no cap. Since the `Graph` constructor rejects
self-loops, such input is never a valid graph and `AcyclicGraphQ` returns `False`; any non-graph
argument also returns `False` (it is a predicate), never unevaluated.
