---
source: src/graph/graphprops.c
---
**Algorithm.** `builtin_empty_graph_q` is a direct structural predicate: it returns `True` iff
its argument is a valid graph whose edge list is empty. The vertex count is irrelevant, so a
graph of isolated vertices counts as empty. No adjacency is built; it reads the length of the
edge `List` (the graph's second argument) after validation.

**Data structures.** None beyond the validated-graph node — the edge-list length comes straight
from the expression, and validation itself is served from the per-node memo.

**Complexity / limits.** `O(1)` after the graph's one-time validation. Exactly one argument. A
non-graph argument returns `False` (it is a predicate), never unevaluated.
