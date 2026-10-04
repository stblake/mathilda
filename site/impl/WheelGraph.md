---
source: src/graph/gmet_generators.c
---
**Algorithm.** `builtin_wheel_graph[n]` joins hub vertex 1 to every other vertex (`n - 1` spokes) and, for `n >= 4`, adds the rim cycle on vertices `2..n`: edges `2-3`, `3-4`, ..., `(n-1)-n`, and the closing edge `2-n`. `WheelGraph[1]` is a single vertex. `n = 2` and `n = 3` would be multigraphs and are left unevaluated, as in Mathematica.

**Data structures.** Spokes and rim edges are appended to the packed-64-bit `Pairs` buffer in already-sorted order, then `pairs_graph` builds `Graph[Range[n], {UndirectedEdge[i, j], ...}]` with a shared `UndirectedEdge` head node and shared vertex integers.

**Complexity / limits.** `O(n)` time, `2(n - 1)` edges for `n >= 4`. Non-integer, non-positive, or over 10^8 vertex arguments leave the call unevaluated; options are not supported.
