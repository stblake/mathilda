---
source: src/graph/directedq.c
---
**Algorithm.** `builtin_directed_graph_q` is `True` exactly when the graph is
valid, has at least one edge, and *every* edge is a `DirectedEdge`. It reads the
directed-edge count from the validated-graph memo (`graph_directed_edge_count`,
`O(1)` after the first query on the node) and compares it with the total edge
count: equal and nonzero means fully directed. An edgeless graph counts as
undirected (as in the Wolfram Language), so `DirectedGraphQ` and
`UndirectedGraphQ` are never both `True`; a mixed graph fails both.

**Data structures.** None of its own — the whole test is two integers read off
the memo entry that `graph_is_valid` populates (`src/graph/graph_util.c`). The
result is a fresh `True`/`False` symbol; the evaluator frees the argument.

**Complexity / limits.** `O(1)` given a memoized graph (`O(V + E)` to populate
the memo on first contact). A non-graph argument gives `False`, never
unevaluated.
