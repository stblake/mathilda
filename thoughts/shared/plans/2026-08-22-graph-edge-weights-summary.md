---
created: 2026-08-22
type: plan
lifecycle: active
status: approved
full_plan: thoughts/shared/plans/2026-08-22-graph-edge-weights.md
---

# Plan Summary: Graph Edge Weights + WeightedAdjacencyMatrix

**Full plan (appendix)**: `thoughts/shared/plans/2026-08-22-graph-edge-weights.md`

## Recommendation
Add `Graph[v, e, EdgeWeight -> {...}]`, `EdgeWeight[g]`, and `WeightedAdjacencyMatrix[g]` —
the one extension `src/graph/adjmat.c:9-10` already flags as pre-approved future work.
Additive only; existing 2-arg `Graph[v,e]` behavior is unchanged.

## Options Considered
1. **3rd constructor arg, `EdgeWeight -> List[...]`** (chosen) — additive, doesn't reuse the
   already-rejected 3-argument-edge shape, matches the codebase's existing trailing-`Rule`
   option idiom.
2. **3-argument edge head** (`DirectedEdge[u,v,w]`) — rejected: already documented and
   enforced as malformed input.
3. **Global `Options`/`OptionValue` registration** — rejected: designed for symbols with
   independently queryable/settable defaults, not a one-shot constructor argument.

## Decisions
- Weights match edges by position in the given edge list.
- `EdgeWeight[g]` defaults to all-`1`s when unweighted.
- No packed/NDArray/`Compile[]` support — documented exemption, consistent with existing
  `AdjacencyMatrix`/`IncidenceMatrix` precedent.
- No weighted `FindShortestPath`/`GraphDistance` in this pass (deferred).
- No derived-vertex weighted construction (`Graph[e, EdgeWeight->w]`) — explicit-vertex form
  only (added after plan-reviewer flagged the silent scope gap).

## Non-goals
Vertex weights, multigraphs, hypergraphs, edge tags beyond `EdgeWeight`, `GraphPlot`
rendering changes, packed/NDArray/`Compile[]` support, weighted shortest-path, derived-vertex
weighted construction.

## Open Questions

### Unresolved
_None._

### Resolved
See full plan — all resolved during research and the plan-reviewer pass.

## Requires Approval
_None._ — scope confirmed directly with the maintainer during research.

## Architecture Impact
- New services introduced: none
- APIs changed: none
- Data crossing a service boundary: none
- New external dependency: none
- Deployment topology change: none

## Subsystems & Dependencies
- Subsystems touched: graph (invocation: inline)
- Interdependencies surfaced: none

## Plan Review outcome
One BLOCKING finding from the adversarial `plan-reviewer` pass — the plan originally claimed
`graph_is_valid` was the sole validation choke point, but `graph_build_adj`
(`src/graph/graph_util.c:206-209`) is a second, independent one that 8 of the 27 builtins
route through, bypassing `graph_is_valid` entirely. Verified directly against source, then
fixed: both choke points now widen via a shared helper, with a new AC (AC-11) and manual
verification step covering all 8 affected builtins. Two WORTH FLAGGING findings (a wrong
line citation, an unstated scope boundary on derived-vertex weighted construction) also
fixed. See the full plan's `## Plan Review` section for the complete transcript.
