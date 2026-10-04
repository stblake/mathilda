---
references:
  - "U. Brandes, *A faster algorithm for betweenness centrality*, J. Math. Sociology **25** (2001) 163-177."
source: src/graph/gmet_centrality.c
---
**Algorithm.** `builtin_betweenness_centrality` computes the (unnormalised) Freeman
betweenness of every vertex with **Brandes' algorithm**: one shortest-path exploration per
source accumulates each vertex's dependency `delta` by walking the settled vertices back in
reverse order of discovery, so the pair-counting is done without materialising predecessor
lists. The vertex form deliberately ignores `EdgeWeight` — it builds its CSR with
`NULL` weights and always takes the unweighted BFS variant (`brandes_bfs_work`); the sibling
edge form accumulates onto edge ids instead and does use the Dijkstra variant
(`brandes_dijkstra_work`, binary heap `GmetHeap`) when the graph carries usable weights. A
mixed directed/undirected graph is left unevaluated, because no standard betweenness
definition fits it.

**Data structures.** Its own compressed-row adjacency `GmetCSR` (`off[]`/`adj[]`, out-arcs),
not the shared `GraphAdj`. Per source it fills `dist[]` (BFS layers), `sigma[]` (shortest-path
counts, as `double`), `delta[]` (dependency), and an `order[]` stack of settled vertices.
Sources are split across a thread team (up to 32) with per-thread accumulators merged at the
end. An undirected graph counts each unordered pair once (a final ×0.5); a directed graph
counts ordered pairs. Results come back as a packed machine-real vector.

**Complexity / limits.** `O(nm)` for the unweighted (BFS) path and `O(nm + n^2 log n)` for the
weighted Dijkstra path, `n = |V|`, `m = |E|`. No size cap. Non-graph or mixed input returns
unevaluated.
