---
source: src/graph/gmet_centrality.c
---
**Algorithm.** `builtin_closeness_centrality` reduces from a cached per-source distance
summary, `gmet_dist_summary(g)`, which runs a multi-source traversal — BFS when the graph is
unweighted, Dijkstra when it carries usable `EdgeWeight` lengths — from every vertex. For each
vertex `v` it records `reach[v]`, the number of vertices `v` can reach, and `sum[v]`, the total
distance to those vertices. The centrality is then `v = reach[v] / sum[v]` (and `0` when `v`
reaches nothing). Dividing by the sum over only the *reachable* set keeps the value finite on a
disconnected graph, which is the practical difference from the textbook `(n-1)/sum` definition.

**Data structures.** The shared summary holds the `reach[]` and `sum[]` arrays (and the
eccentricity used by `GraphCenter`/`GraphRadius`), memoised on the graph node so repeated metric
queries on the same graph are `O(1)`. Distances follow edge direction on a directed graph and an
undirected edge is usable both ways. The output is a packed machine-real vector.

**Complexity / limits.** The summary costs `O(V(V+E))` unweighted and `O(V(E + V log V))`
weighted; the reduction to the centrality vector is `O(V)`. Single argument only; non-graph
input, or a summary that fails to build, returns unevaluated.
