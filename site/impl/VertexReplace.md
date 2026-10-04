---
source: src/graph/gops_edit.c
---
**Algorithm.** `builtin_vertex_replace` takes `VertexReplace[g, rule]` or `VertexReplace[g, {rules}]`, where every rule must be a `Rule` or `RuleDelayed`. It evaluates `Replace[VertexList[g], rules, {1}]` once, so first-match-wins, patterns and delayed right-hand sides behave exactly as in `Replace`. The images become the new vertex list, de-duplicated through a scratch hash with a position map from old to new index. Each edge is rebuilt with its original head and orientation from the mapped endpoints (untouched edges are shared), and weights are copied.

**Data structures.** `Graph[List, List]` expression tree; the old-to-new vertex map is an `int` array, and the rebuilt endpoint arrays go to `gops_graph_new`.

**Complexity / limits.** `O(V + E)` plus one `Replace` evaluation over the vertex list. If the replacement merges two vertices so that a self-loop or parallel edge would result (for example renaming 2 to 1 on a path 1-2-3), the graph is not simple and the call stays unevaluated. Non-rule second arguments are also left alone.
