---
created: 2026-08-23T04:15:30Z
researcher: Michael Sollami
topic: "Second graph extension: weighted shortest path"
type: research
lifecycle: active
full_research: thoughts/shared/research/2026-08-23-weighted-shortest-path.md
---

# Research Summary: Weighted Shortest Path

**Full research (appendix)**: `thoughts/shared/research/2026-08-23-weighted-shortest-path.md`

## Recommendation
Make `FindShortestPath`/`GraphDistance` weight-aware (Dijkstra) when the graph carries
non-negative numeric `EdgeWeight`, matching real Wolfram Language semantics and closing the
exact gap ticket 1's own Non-goals named as deferred follow-up work. Build a local,
call-scoped weighted adjacency rather than extending the shared `GraphAdj` structure 8 other
builtins depend on.

## Options Considered
1. **Dijkstra dispatch inside the existing two builtins, local weighted structure** (chosen)
   — matches real WL behavior, minimal blast radius, reuses ticket 1's
   `graph_resolve_edge_weights`.
2. **Extend `GraphAdj` itself with a weight array** — rejected: widens a structure 8 builtins
   share, the exact shape of risk ticket 1's `plan-reviewer` pass flagged as a real defect.
3. **New, separately-named builtins (`WeightedFindShortestPath`)** — rejected: not how real
   Wolfram Language behaves (weight-awareness is automatic based on graph properties, not a
   separate function name), and duplicates two builtins instead of completing them.

## Decisions
- Dijkstra fires only when `EdgeWeight` is present AND every weight is non-negative numeric;
  otherwise, falls back to the existing unweighted BFS exactly as before (fails safe, never
  regresses a previously-working call to unevaluated).
- No change to `GraphAdj`/`graph_build_adj` — a separate, local structure only where needed.

## Non-goals
No Bellman-Ford / negative-weight support. No change to any other `graph_build_adj`-routed
builtin. No A*/bidirectional search — plain O(V²) Dijkstra, matching this codebase's own
complexity tolerance for small-graph exact algorithms (`VertexConnectivity`'s precedent).

## Open Questions

### Unresolved
_None._

### Resolved
See full research — all three resolved directly from the codebase, no maintainer
consultation needed for this pass beyond what ticket 1 already decided.

## Requires Approval
None — this is the explicit, named follow-up from ticket 1's own Non-goals.
