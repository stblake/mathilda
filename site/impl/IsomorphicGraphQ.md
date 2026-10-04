---
references:
  - "B. D. McKay and A. Piperno, *Practical graph isomorphism, II*, J. Symbolic Comput. **60** (2014) 94-112."
  - "R. Paige and R. E. Tarjan, *Three partition refinement algorithms*, SIAM J. Comput. **16** (1987) 973-989."
source: src/graph/galg_iso.c
---
**Algorithm.** `builtin_isomorphic_graph_q` decides whether all its arguments are pairwise
isomorphic by an individualization-refinement canonical-labelling engine, the family of nauty /
bliss / Traces. The core is colour refinement — one-dimensional Weisfeiler–Leman implemented in
the Hopcroft / Paige–Tarjan style (a FIFO of splitter cells, skipping the largest fragment for
an `O((n+m) log n)` pass). A running 64-bit "trace" hash is emitted per split and compared
on-the-fly against a reference, so a branch aborts the moment it can no longer match. Over that
sits an individualization-refinement search tree with exact per-level undo, automorphism
detection and orbit pruning. The two graphs are compared first by a lockstep `G -> H` trace
search; if that exceeds a node limit it falls back to comparing the two canonical forms. Cheap
rejections (vertex/edge counts, colour-class sizes, degree histograms) run first.

**Data structures.** `GalgIsoGraph` holds up to three CSR relations — undirected, out and in —
plus optional vertex colours, so directed graphs, loops and multigraphs are handled via a
coloured subdivision. The ordered partition carries label/inverse arrays, cell start/size
arrays, and a min segment tree over positions for target-cell selection; automorphisms are kept
as sparse generators with union-find orbits.

**Complexity / limits.** One refinement is `O((n+m) log n)`; the search tree is small for random
or structured graphs, with a worst-case refinement-node budget `GI_BUDGET = 5·10^7` (polled
against `TimeConstrained`). Variadic, at least two arguments; `True` iff all are mutually
isomorphic. Any non-graph argument gives `False`; only budget exhaustion or allocation failure
returns unevaluated.
