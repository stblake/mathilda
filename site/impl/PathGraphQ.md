---
source: src/graph/gops_preds.c
---
**Algorithm.** `builtin_path_graph_q` follows Mathematica's definition of a path:
at least one vertex, connected, and with bounded degrees — undirected, every
degree `<= 2`; directed, every in- and out-degree `<= 1`. One pass over the edges
accumulates in/out degrees and fails as soon as a bound is exceeded; a union-find
over the edges then checks connectivity (the number of joins must equal `n - 1`).
Because a cycle meets the degree bounds and is connected, it counts as a (closed)
path — `PathGraphQ[CycleGraph[3]]` is `True`, as in Mathematica 15. A mixed graph
is never a path, and a graph with more than `n` edges fails the degree bounds
immediately.

**Data structures.** A `GopsView` plus three `O(V)` arrays: in-degree `din[]`,
out-degree `dout[]`, and the union-find `parent[]`. All work is on the view's
integer endpoint arrays.

**Complexity / limits.** `O(V + E)`. A non-graph argument gives `False`. The
early `ne > n` rejection keeps a dense graph from even allocating the degree
arrays.
