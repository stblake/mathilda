---
created: 2026-08-23T04:15:30Z
researcher: Michael Sollami
source_sha: 3d87224771270a91d7402de3f8b9908dd33e58c8
branch: main
repository: mathilda
topic: "Second graph extension: what does the code/tests show a maintainer would want next?"
tags: [research, codebase, graph, shortest-path, dijkstra]
subsystems: [graph]
type: research
lifecycle: active
status: complete
last_updated: 2026-08-23
last_updated_by: Michael Sollami
---

# Research: Second graph extension — weighted shortest path

**Date**: 2026-08-23T04:15:30Z
**Researcher**: Michael Sollami
**Git Commit**: 3d87224771270a91d7402de3f8b9908dd33e58c8
**Branch**: main
**Repository**: mathilda

## TL;DR
Ticket 1 (edge weights) explicitly deferred making `FindShortestPath`/`GraphDistance`
weight-aware, naming it in its own `Non-goals` as the natural next step. `GraphAdj`
(`graph_util.c`, shared by 8 builtins) stores **no per-edge weight at all** — only
successor/predecessor vertex indices — so real Wolfram Language semantics (both builtins
auto-dispatch to a weighted algorithm when `EdgeWeight` is present) require a real Dijkstra
implementation, not a config flag. This is harder than ticket 1: a new algorithm, not just a
new builtin, and it changes two existing builtins' behavior on weighted graphs rather than
adding new read-only ones. Sized at a few hours given a simple O(V²) array-based Dijkstra
(matching this codebase's own precedent: `VertexConnectivity`'s docstring calls itself
"exact brute-force ... intended for small graphs").

## Summary
`src/graph/shortestpath.c` implements unweighted BFS for both `FindShortestPath[g,s,t]` and
`GraphDistance[g,s,t]`, routed through the shared `GraphAdj` (`graph_build_adj`). Real
Wolfram Language's own `FindShortestPath`/`GraphDistance` automatically use edge weights when
present and fall back to unweighted BFS otherwise — that is the behavior ticket 1's own
research and plan explicitly named as deferred (`thoughts/shared/research/
2026-08-22-graph-edge-weights-extension.md`'s Non-goals: *"No weighted-shortest-path /
Dijkstra mode ... deferred to a follow-up"*). `GraphAdj` has no weight storage, so this
requires building a small, local weighted-adjacency pass (reusing `graph_resolve_edge_weights`
from ticket 1) rather than touching the shared structure 8 other builtins depend on —
learning directly from ticket 1's plan-reviewer-caught lesson about `graph_build_adj` being a
sensitive shared choke point.

## Open Questions

### Unresolved
_None._

### Resolved
- [x] Does `GraphAdj` already carry weights that a Dijkstra pass could reuse? — No, confirmed
      by reading `graph_util.c`'s `GraphAdj` struct and `graph_build_adj`'s fill loop: only
      `int` successor/predecessor indices, no weight field anywhere.
  - [x] Should this touch the shared `GraphAdj`/`graph_build_adj`? — No: build a local,
      call-scoped weighted adjacency inside `shortestpath.c` instead, to avoid widening the
      blast radius of a structure 8 builtins depend on (direct lesson from ticket 1's
      `plan-reviewer` finding).
- [x] How should non-numeric or negative weights be handled? — Fall back to the existing
      unweighted BFS behavior rather than failing: Dijkstra requires non-negative numeric
      weights to be correct, and a previously-working call should not start returning
      unevaluated just because a graph happens to carry a symbolic or negative weight.
      Documented as an explicit limitation, matching this codebase's existing style
      (`VertexConnectivity`'s own "intended for small graphs" self-limitation).
- [x] Prior attempt or known constraint? — None found in git history (`git log --oneline --
      src/graph/shortestpath.c` shows only the original subsystem-add commit); this session's
      own ticket-1 Non-goals is the only prior signal, and it points at doing exactly this.

## Requires Approval
None beyond what ticket 1 already flagged and deferred to "a follow-up" — this is that
follow-up.

---

## Research Question
"A second, harder ticket in this repo, chosen the same way as the first — from what the
code and tests show a maintainer would want, sized at a few hours."

## Detailed Findings

### Current implementation (`src/graph/shortestpath.c`, full file read)
- `bfs()` (lines 20-33): unweighted BFS over `GraphAdj.out[]`, O(V+E).
- `resolve()` (lines 36-45): validates `g`, builds `GraphAdj` via `graph_build_adj`,
  resolves `s`/`t` to vertex indices via `graph_vertex_index`.
- `builtin_find_shortest_path`/`builtin_graph_distance`: call `resolve` + `bfs`, no weight
  awareness anywhere.

### `GraphAdj` has no weight field
`graph_util.c`'s `GraphAdj` struct (declared in `graph.h:112-117`): `n`, `verts`
(borrowed), `outdeg`/`out`, `indeg`/`in` — all vertex-index arrays, no weight storage. Adding
a weight array here would touch `graph_build_adj`, used directly by `components.c` (x2),
`connectivity.c` (x2), `spanningtree.c`, and `shortestpath.c` (x2) — 7 call sites across 5
files. Given ticket 1's `plan-reviewer` finding was specifically about this function being a
shared, easy-to-miss choke point, the lower-risk design keeps `graph_build_adj` completely
unchanged and builds a separate, call-scoped weighted structure only where Dijkstra is
actually used.

### Weight resolution is already solved (ticket 1)
`graph_resolve_edge_weights(g)` (`graph_util.c`, added in ticket 1) already returns the
per-edge weight list in `EdgeList` order, defaulting to all-`1`s when unweighted — exactly
the input a Dijkstra pass needs, keyed against `EdgeList[g]`/`edges` the same way
`wtadjmat.c` already consumes it.

### Test coverage
`tests/test_graph.c`'s `test_shortest_path` (added in the original subsystem commit) and
`test_edge_weights` (added in ticket 1, this session) both exist. `test_edge_weights`'s
AC-11 case explicitly asserts `FindShortestPath`/`GraphDistance` **ignore** weights on a
graph carrying them — a real, checked-in test of the exact limitation ticket 1 deliberately
left in place and named as deferred. Implementing this ticket requires **changing that
specific assertion** (not removing the AC-11 pattern — the other 6 `graph_build_adj`-routed
builtins in that test still assert unchanged, still-correct behavior).

## Code References
- `src/graph/shortestpath.c` — the whole file; both builtins and their shared `resolve`/`bfs`
- `src/graph/graph.h:112-120` — `GraphAdj` struct and `graph_build_adj` declaration
- `src/graph/graph_util.c` — `graph_resolve_edge_weights` (ticket 1), reused here
- `tests/test_graph.c`'s `test_edge_weights` — the AC-11 assertions that must change
- `thoughts/shared/research/2026-08-22-graph-edge-weights-extension.md` — ticket 1's
  Non-goals naming this exact follow-up

## Architecture Insights
Same one-builtin-per-file convention; this ticket modifies two existing builtin files rather
than adding new ones, since the change is "make these two smarter," not "add two new
readers." The GraphVIdx O(1) index helper (ticket-1-adjacent infra) is reusable here too.

## Historical Context (from thoughts/)
- `thoughts/shared/research/2026-08-22-graph-edge-weights-extension.md` and
  `thoughts/shared/plans/2026-08-22-graph-edge-weights.md` — ticket 1, whose Non-goals and
  GR-09 plan-reviewer finding both directly shape this ticket's scope and design.

## Related Research
- `thoughts/shared/research/2026-08-22-graph-edge-weights-extension.md`
