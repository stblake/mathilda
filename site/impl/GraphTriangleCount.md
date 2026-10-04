---
references:
  - "N. Chiba and T. Nishizeki, *Arboricity and subgraph listing algorithms*, SIAM J. Comput. **14** (1985) 210-223."
  - "M. Latapy, *Main-memory triangle computations for very large (sparse (power-law)) graphs*, Theoret. Comput. Sci. **407** (2008) 458-473."
source: src/graph/gmet_cluster.c
---
**Algorithm.** `builtin_graph_triangle_count` gives the number of triangles of an
undirected graph, or of directed 3-cycles `u->v->w->u` of a directed graph. The
core `tri_compute` lists triangles over the underlying simple graph with the
**degree-ordered orientation**: each edge is oriented from the endpoint of lower
`(degree, index)` rank to the higher, which bounds every oriented out-degree by
`O(sqrt m)` and the total listing work by `O(m^1.5)` (Chiba-Nishizeki / Latapy).
For each oriented edge `u -> v` it intersects `u`'s and `v`'s forward
neighbourhoods via a mark array. For a directed graph each oriented edge carries
a two-bit record of which of the two arcs exist, so the 3-cycle test is two
bit-ANDs per candidate; a transitive triple counts `0` and a doubly-linked
triangle counts `2`.

**Data structures.** CSR oriented adjacency built from the memo's integer
endpoints: `ooff`/`oadj`, plus a per-edge flag byte `ofl` for the directed case.
Rows are processed on the thread team with per-thread mark arrays and counters,
and the undirected build uses index-order orientation when `maxdeg^2 <= 4m`
(where the two bounds agree) to avoid the degree sort.

**Complexity / limits.** `O(m^1.5)` time, exact `Integer` result. `EdgeWeight` is
ignored; a mixed graph is left unevaluated.
