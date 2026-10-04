---
source: src/graph/graphq.c
---
**Algorithm.** `builtin_graph_q` is a thin wrapper over `graph_is_valid`. It returns `True` for a canonical, valid `Graph[List[vertices], List[edges]]` and `False` for any other expression, including atoms. A call with other than one argument is left unevaluated. Validity means that the edge endpoints are vertices of the graph, that each edge is a `DirectedEdge` or `UndirectedEdge`, and that there are no duplicate or self-loop edges.

**Data structures.** Validation lives in `graph_util.c`, which memoizes the verdict on the graph node together with a vertex hash index, the edge-key set, and the `eu`/`ev`/`directed` edge-index views. Every other graph predicate and algorithm reuses that entry, so `GraphQ` is the call that fills it.

**Complexity / limits.** The first call on a graph is `O(V + E)`. Repeats on the same graph node are an `O(1)` memo hit. It is purely structural and never evaluates or rewrites its argument. Hypergraphs and other graph-like objects are not `Graph` expressions, so `GraphQ` gives `False` for them.
