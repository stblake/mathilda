---
created: 2026-08-23T02:49:49Z
researcher: Michael Sollami
source_sha: 6a9a6259768a0d5516bba7747220d70835a3a79f
branch: main
repository: mathilda
topic: "What genuine, few-hours-sized extension to Graph/HyperGraph would a maintainer want?"
tags: [research, codebase, graph, hypergraph, weighted-graph]
subsystems: [graph]
type: research
lifecycle: active
status: complete
last_updated: 2026-08-23
last_updated_by: Michael Sollami
---

# Research: What genuine, few-hours-sized extension to Graph/HyperGraph would a maintainer want?

**Date**: 2026-08-23T02:49:49Z
**Researcher**: Michael Sollami
**Git Commit**: 6a9a6259768a0d5516bba7747220d70835a3a79f
**Branch**: main
**Repository**: mathilda

## TL;DR
HyperGraph does not exist anywhere in this codebase and is explicitly locked out of MVP
scope alongside multigraphs/parallel-edges/edge-tags — building it properly is a multi-day
effort, not a few hours, and would contradict a documented scope decision. Edge weights and
`WeightedAdjacencyMatrix`, by contrast, are the one thing the source itself flags as a
pre-approved "future hook" (`src/graph/adjmat.c:9-10`, `docs/spec/builtins/graphs.md:19-21`)
and are unimplemented, untested, and unconsumed downstream — the right-sized, justified
choice. Open question: whether weighted variants of `FindShortestPath`/`GraphDistance`
belong in this pass (recommend: no, out of scope).

## Summary
`src/graph/` (19 files, ~1850 LoC, 27 builtins, all `ATTR_PROTECTED`) implements a
Mathematica-style simple-graph subsystem: construction/validation, query/representation,
matrix views, generators, BFS/Tarjan-based search, and a circle-layout `GraphPlot`. Three
commits total (`56035303` initial add, `de2ef6e4` unrelated evaluator perf work,
`8d71d845` a vertex-index perf fix) — no weighted/multigraph/hypergraph exploration has
happened. `docs/spec/builtins/graphs.md:19-21` locks the MVP to simple graphs with **no
hypergraphs, no multigraphs, no edge tags, no edge/vertex weights** — except it explicitly
carves out `WeightedAdjacencyMatrix`/edge weights as a **documented future extension**, and
`src/graph/adjmat.c:9-10`'s own comment says so again in-code. `AdjacencyMatrix`/
`IncidenceMatrix` output is consumed nowhere internally (`src/linalg/`, `src/compile/` never
reference them) — a weighted variant needs no downstream migration. Neither is on
`src/pack.c`'s `AWARE` list; graph builtins are structural (operate on graph trees, not
numeric buffers), so this is consistent precedent for a new matrix-returning builtin to
follow, not a gap this change needs to close.

## Open Questions

### Unresolved
- [ ] Should `FindShortestPath`/`GraphDistance` gain a weighted (Dijkstra) mode in the same
      pass, now that weights exist? Recommend deferring — real scope growth (a new
      algorithm, not a new builtin) beyond "a few hours."

### Resolved
- [x] Does HyperGraph exist anywhere in this codebase? — No (exhaustive grep across
      `src/`, `tests/`, `docs/`, `src/internal/*.m`); it is named only as an explicitly
      out-of-scope MVP exclusion in `docs/spec/builtins/graphs.md:20`.
- [x] Is there a pre-approved "next" extension already flagged in the code? — Yes:
      `WeightedAdjacencyMatrix`/edge weights, `src/graph/adjmat.c:9-10` and
      `docs/spec/builtins/graphs.md:19-21`.
- [x] Do any consumers assume `AdjacencyMatrix`/`IncidenceMatrix` are exactly 0/1/-1? — No;
      grepped all of `src/` and `src/internal/*.m` — no internal consumer exists at all.
- [x] Does the packed/NDArray/Compile[] mandate (CLAUDE.md) apply here? — No current graph
      builtin (including the existing `AdjacencyMatrix`) is on `src/pack.c`'s `AWARE` list;
      they are structural/constructive, not elementwise-numeric. A new matrix-returning
      builtin follows the same, already-established precedent.
- [x] Any prior attempt or deliberate constraint the maintainer already knows about? —
      Asked directly (grill-me, research-open); answer: no, research fresh from the code.

## Requires Approval
Whether to scope this to edge weights + `WeightedAdjacencyMatrix` + an `EdgeWeight[g]`
query builtin only (recommended), versus also touching `FindShortestPath`/`GraphDistance`
to make them weight-aware. Recommend the narrower scope; flagged above as Unresolved.

---

## Research Question
"Extend the Graph and HyperGraph functionality in this codebase. Pick a genuine, useful
extension that a maintainer would actually want, sized at a few hours rather than a week,
justified from what the code and tests actually show."

## Detailed Findings

### Graph subsystem inventory
- `src/graph/graph.c` registers 27 builtins (`Graph`, `GraphQ`, `VertexList`, `EdgeList`,
  `VertexCount`, `EdgeCount`, `AdjacencyList`, `VertexDegree`, `VertexInDegree`,
  `VertexOutDegree`, `DirectedGraphQ`, `AdjacencyMatrix`, `IncidenceMatrix`,
  `AdjacencyGraph`, `CompleteGraph`, `CycleGraph`, `PathGraph`, `RandomGraph`,
  `FindShortestPath`, `GraphDistance`, `ConnectedComponents`,
  `WeaklyConnectedComponents`, `StronglyConnectedComponents`, `FindSpanningTree`,
  `ConnectedGraphQ`, `VertexConnectivity`, `GraphPlot`), each `symtab_add_builtin` call
  paired with `symtab_set_docstring` and `ATTR_PROTECTED`.
- One builtin per `.c` file (mirrors `src/linalg/`'s layout), per `src/graph/graph.h:1-22`'s
  header comment.
- Canonical form: `Graph[List[verts], List[edges]]`, edges `DirectedEdge[u,v]` /
  `UndirectedEdge[u,v]`, `Rule`/`TwoWayRule` accepted as construction-time sugar
  (`src/graph/graph.h:8-13`).

### Locked scope and the one documented exception
`docs/spec/builtins/graphs.md:19-21`:
> "MVP scope (locked): simple graphs only — no parallel edges, no self-loops, no edge tags,
> no multigraphs, no hypergraphs, and no edge/vertex weights. `WeightedAdjacencyMatrix` and
> edge weights are a documented future extension."

`src/graph/adjmat.c:9-10` (in-code, right next to the current `AdjacencyMatrix`
implementation):
> "Future hook: a WeightedAdjacencyMatrix would fill entries with edge weights instead of 1
> (Locked Decision 2); not implemented in the MVP."

No other TODO/FIXME/"not implemented" comment exists anywhere in `src/graph/*.c`. The only
other named future hooks are cosmetic: arrowheads for directed edges in `GraphPlot`
(`docs/spec/builtins/graphs.md:156-157`) and a force-directed layout
(`src/graph/graphplot.c`), neither of which changes graph semantics or query results.

### Test coverage (`tests/test_graph.c`, 319 lines, 15 test functions, all passing at
`source_sha` above — verified: `./tests/build-main/graph_tests` → "All graph tests passed!")
Every one of the 27 builtins has at least one asserting test. Thinnest coverage:
`VertexConnectivity` (4 hardcoded cases: `PathGraph[4]→1`, `CycleGraph[5]→2`,
`CompleteGraph[4]→3`, a disconnected pair→0) and boundary shapes (single-vertex,
zero-edge, disconnected graphs) get one or two cases each rather than systematic coverage.
No existing test constructs a graph with a third, option-style constructor argument — the
constructor's argument-count handling is exactly `Graph[v,e]`/`Graph[e]` today
(`src/graph/construct.c`), which is the surface a weighted-edge extension has to grow
without breaking.

### Downstream consumption of matrix output
`AdjacencyMatrix`/`IncidenceMatrix` are grepped across the entire `src/` tree and
`src/internal/*.m`: the only hits are the symbol interning table
(`src/sym_names.{h,c}`) and the graph subsystem's own files. No internal consumer in
`src/linalg/`, `src/compile/`, or elsewheer assumes the matrix entries are literally `0`/`1`
(or `-1/0/1` for the incidence form) — a real-valued weighted variant requires zero changes
outside the graph subsystem.

### Packed/NDArray/Compile[] applicability
`src/pack.c`'s `AWARE` list (the CLAUDE.md-mandated registry of heads with a numeric fast
path) contains no `Graph`/`Adjacency`/`Vertex`/`Edge` entry — confirmed by grep — and
neither does its `NOT_AWARE`/`INT64_OK` commentary. This is not a gap introduced by this
research: `AdjacencyMatrix`, already shipped, is in the same position. Graph builtins
consume/produce compound `Expr` trees (graphs, vertex/edge lists), not numeric buffers —
CLAUDE.md's own escape hatch ("a purely symbolic/structural head... genuinely cannot
support a surface") applies to the whole subsystem as shipped. A new `WeightedAdjacencyMatrix`
follows the identical shape as `AdjacencyMatrix` (builds a fresh matrix from graph structure)
and inherits the same, already-established position — not a new exemption to invent.

### Git history
Three commits touch `src/graph/`: `56035303` (initial subsystem, 2026-06-29), `de2ef6e4`
(unrelated evaluator-wide perf work, 2026-07-30), `8d71d845` (vertex-index perf fix,
2026-08-05, motivated by a 20000-vertex benchmark). No weighted-graph or hypergraph work has
been attempted or reverted.

## Code References
- `src/graph/graph.h:1-22` — subsystem header comment, canonical form, ownership contract
- `src/graph/graph.c` — builtin registration, attributes, docstrings (all 27 builtins)
- `src/graph/construct.c` — `Graph[...]` construction/validation/normalization
- `src/graph/adjmat.c:9-10` — the in-code `WeightedAdjacencyMatrix` future-hook comment
- `docs/spec/builtins/graphs.md:19-21` — the locked MVP-scope paragraph
- `tests/test_graph.c` — the full existing test suite (baseline: all 15 tests pass)
- `src/pack.c` — `AWARE`/`NOT_AWARE`/`INT64_OK` registries (no graph entries either way)

## Architecture Insights
- One-builtin-per-file, hub-registers-in-`graph.c` is the fixed convention (mirrors
  `src/linalg/`); a new builtin follows the same shape.
- Options-as-trailing-`Rule[]` is the established idiom for optional builtin arguments
  elsewhere in the codebase (`src/numerical_calculus/nderiv.c:531-566`,
  `src/numerical_calculus/ndsolve.c:90-109` — a local `is_option`/`apply_option` pair per
  module), not the global `Options`/`OptionValue`/`SetOptions` registry
  (`src/options_builtin.c`), which is for symbols whose options are queried/set
  independently of a single call. `Graph[v, e, EdgeWeight -> {...}]` should follow the
  local trailing-rule idiom, matching how the constructor already accepts `Rule`/
  `TwoWayRule` edge sugar positionally.
- Vertex/edge lookups already have an O(1) hash-index helper (`GraphVIdx`,
  `src/graph/graph_util.c`) introduced specifically to keep new graph builtins from
  reintroducing the O(E·V) scans fixed in `8d71d845` — any new builtin should build on it
  rather than re-scanning `List` args linearly.

## Historical Context (from thoughts/)
None — `thoughts/` did not exist before this research (created for `NOTES_DIR` output).

## Related Research
None yet — this is the first research document under `thoughts/shared/research/`.
