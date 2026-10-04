---
references:
  - "C. Hierholzer, *Über die Möglichkeit, einen Linienzug ohne Wiederholung und ohne Unterbrechung zu umfahren*, Mathematische Annalen **6** (1873) 30-32."
source: src/graph/gops_cycles.c
---
**Algorithm.** `builtin_find_eulerian_cycle` returns `{cycle}` — a closed walk
using every edge exactly once, as a list of edges — or `{}` when the graph has no
Eulerian cycle. It first calls `gops_eulerian`, which decides the precondition
from degrees and connectivity alone: for a directed graph every vertex must have
in-degree equal to out-degree, for an undirected graph every vertex must have even
degree, and the edge-touched vertices must form one connected block (balance plus
weak connectivity gives strong connectivity in the directed case). A non-Eulerian
graph gives `{}`; an edgeless graph with at least one vertex gives `{{}}`; the null
graph `{}`. When the graph is Eulerian the walk is built by **Hierholzer's
algorithm**, run iteratively: starting at the first vertex in `VertexList` order
that has an incident edge, it follows unused edges (taken in `EdgeList` order)
until it returns stuck, popping vertices to record the tour. Directed cycles are
reported forwards; undirected cycles in Hierholzer's pop order, each undirected
edge written in the direction it was actually walked, as Mathematica does. Only
`FindEulerianCycle[g]` and `FindEulerianCycle[g, 1]` are supported — mixed graphs,
and a second argument other than `1` (`n > 1` or `All`), leave the call unevaluated.

**Data structures.** The graph is read through a `GopsView`; traversal uses a
cached incidence CSR (`GOPS_INC_OUT` for a directed graph, `GOPS_INC_ALL` for an
undirected one). Scratch is a `used[]` edge-visited bitmap, a per-vertex incidence
cursor `ptr[]`, an explicit DFS stack (`sv`/`se` for vertex and entering edge), and
pop-order arrays (`cv`/`ce`) — no recursion, so there is no call-stack bound. The
result is `{List of edges}`, each edge `expr_copy`'d from the original (or
re-oriented for a reversed undirected traversal via `walked_edge`).

**Complexity / limits.** `O(V + E)` in both the Eulerian test and the walk, with
one incidence cursor advance per edge. A `TimeConstrained` abort can `longjmp` out
and leak the scratch arrays, the standard tradeoff across `src/graph/`. The output
is a concrete edge list, not a graph.
