---
references:
  - "B. D. McKay and A. Piperno, *Practical graph isomorphism, II*, J. Symbolic Comput. **60** (2014) 94-112."
  - "B. Weisfeiler and A. Leman, *The reduction of a graph to canonical form and the algebra which appears therein*, Nauchno-Techn. Inform. Ser. 2, **9** (1968) 12-16."
source: src/graph/galg_isoheads.c
---
**Algorithm.** `builtin_canonical_graph` returns a relabelling of `g` onto the
vertices `1..n` such that two graphs are isomorphic iff their canonical graphs
are identical. It first reduces `g` to a vertex-coloured simple structure
(`gi_build`): self-loops fold into the vertex colour, and an edge class of
multiplicity `k > 1` is subdivided through a new colour-`(k, kind)` vertex, so
mixed graphs, loops and multigraphs are all handled. `gi_canonical_rank` then
runs `galg_iso_canon` — the individualization-refinement engine of
`galg_iso.c`, the nauty / bliss / Traces family: a 1-dimensional Weisfeiler-Leman
colour refinement (Hopcroft / Paige-Tarjan splitter queue) to the coarsest
equitable partition, a search tree that individualizes one vertex of the target
cell per level, and a leaf key `(trace₁, …, relabelled graph)` whose best leaf
defines the canonical labelling. The ranks of the original vertices are read off
the canonical leaf, every edge is rewritten on `1..n` (undirected endpoints
sorted), the edge list is sorted, and a fresh `Graph[{1,…,n}, {…}]` is built.

**Data structures.** `GalgIsoGraph` carries up to three CSR relations
(undirected, directed-out, directed-in) plus a ranked `vcol` colour array; the
engine is pure C with no `Expr` dependency. The ordered partition is kept
nauty-style in a `lab[]` position array (a cell named by its start), refinement
is driven by a FIFO of splitter cells, and each refinement emits a 64-bit
**trace** compared on the fly so a losing branch is abandoned early. Edges are
sorted with `qsort` over `(a, b, directed)` triples.

**Complexity / limits.** One refinement is `O((n+m) log n)`; the search tree is
small for almost all graphs but can be super-polynomial on adversarial
constructions (CFI-type), so the whole search is budgeted at `GI_BUDGET = 5·10⁷`
nodes and polls the `TimeConstrained` deadline — on exhaustion the head stays
unevaluated rather than guessing. Edge weights and other properties are ignored,
as in Mathematica.
