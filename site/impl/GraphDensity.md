---
source: src/graph/gmet_distance.c
---
**Algorithm.** `builtin_graph_density` reads three numbers straight off the `Graph`
expression and does no traversal. It counts the directed arcs `nd` with
`graph_directed_edge_count`, reads the total edge count `ne` and the vertex count
`n` from the lengths of the graph's vertex-list and edge-list arguments, and
returns the exact rational `(nd + 2(ne − nd)) / (n(n−1))`. Each directed arc
contributes one to the numerator and each undirected edge two, so the numerator is
the number of ordered endpoint pairs joined by an edge; the denominator
`n(n−1)` is the number of ordered pairs of distinct vertices, i.e. the maximum
possible directed arcs on `n` labelled vertices. A complete undirected graph
therefore has density `1`.

**Data structures.** None beyond the input `Expr` tree: the arg counts are read
directly, and `make_rational` builds the result, reducing it to lowest terms.

**Complexity / limits.** `O(E)` to classify the edges as directed or undirected,
then `O(1)` arithmetic; the answer is an exact `Rational` with no floating point.
Graphs with fewer than two vertices leave the call unevaluated (the density is
undefined — `n(n−1) = 0`), matching Wolfram. Self-loops, if present, are counted
among the edges by `ne`.
