---
source: src/graph/gops_setops.c
---
**Algorithm.** `builtin_graph_intersection` builds the graph on the union of the argument graphs'
vertex sets whose edges are those common to *every* graph. Vertices are merged in canonical
(`expr_compare`) order. The edges of the first graph are tested for membership in each other
graph — an edge is kept when it is present in all of them (`hits == ng-1`) — and the kept edges
are counting-sorted into canonical order. The first graph's orientation is preserved, and edge
weights are dropped. `GraphIntersection[g]` of a single graph returns `g`.

**Data structures.** A `GopsView` per graph with integer endpoints; a `GopsKeySet` hash of the
first graph's edges, probed by each other graph's edges mapped onto the shared vertex indices;
the vertex union sort uses radix/counting for all-integer vertices and a comparison sort
otherwise.

**Complexity / limits.** `O(V + E)` apart from the `O(V log V)` vertex-union sort (`O(V)` when
vertices are already ordered integers). Variadic, at least one argument; a non-graph argument
returns unevaluated. An undirected edge equals its reversal under the edge key, so orientation
does not block a match.
