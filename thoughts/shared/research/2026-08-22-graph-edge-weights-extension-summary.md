---
created: 2026-08-23T02:49:49Z
researcher: Michael Sollami
topic: "What genuine, few-hours-sized extension to Graph/HyperGraph would a maintainer want?"
type: research
lifecycle: active
full_research: thoughts/shared/research/2026-08-22-graph-edge-weights-extension.md
---

# Research Summary: What genuine, few-hours-sized extension to Graph/HyperGraph would a maintainer want?

**Full research (appendix)**: `thoughts/shared/research/2026-08-22-graph-edge-weights-extension.md`

## Recommendation
Add edge weights + `WeightedAdjacencyMatrix` (+ a matching `EdgeWeight[g]` query builtin).
This is the one extension the codebase itself already names as pre-approved future work
(`src/graph/adjmat.c:9-10`, `docs/spec/builtins/graphs.md:19-21`), it's unimplemented,
untested, and has zero downstream consumers to migrate — right-sized for a few hours.
HyperGraph does not exist and is explicitly locked out of MVP scope; building it properly
(construction, validation, printing, its own query/generator/algorithm surface) is a
multi-day effort and would override a documented scope decision, not extend it.

## Options Considered
1. **Edge weights + `WeightedAdjacencyMatrix` + `EdgeWeight[g]`** — pre-approved in-code,
   zero downstream migration, self-contained to `src/graph/`. Chosen.
2. **HyperGraph from scratch** — explicitly locked out of MVP scope
   (`docs/spec/builtins/graphs.md:20`); needs its own construction/validation/printing/
   query/generator surface — a week-plus effort, not a few hours.
3. **`GraphPlot` arrowheads / force-directed layout** — real gaps, but cosmetic/visual only;
   doesn't touch graph semantics or add queryable functionality.
4. **`VertexConnectivity` test hardening** — a legitimate but much smaller test-only task,
   not really "extending functionality."

## Decision Criteria
- Must be justified by what the code/docs/tests actually show, not convenience: option 1 is
  the only one the source explicitly names as intended future work.
- Must fit "a few hours": option 1 is one new small builtin plus a constructor extension
  that's additive to the existing 2-arg canonical form; option 2 fails this outright.
- Must not contradict a locked scope decision: option 2 would.

## Open Questions

### Unresolved
- [ ] Should `FindShortestPath`/`GraphDistance` gain a weighted (Dijkstra) mode in this same
      pass? Recommend deferring to a follow-up — real algorithmic scope growth.

### Resolved
- [x] HyperGraph existence — confirmed absent, confirmed locked out of scope.
- [x] Pre-approved next extension — confirmed: weighted edges / `WeightedAdjacencyMatrix`.
- [x] Downstream consumers of `AdjacencyMatrix`/`IncidenceMatrix` assuming 0/1 values — none.
- [x] Packed/NDArray/Compile[] mandate applicability — none of the existing graph builtins
      are on `src/pack.c`'s `AWARE` list; structural/constructive, consistent precedent.
- [x] Prior attempt / known constraint — none; confirmed directly with the maintainer.

## Requires Approval
Scope: edge weights + `WeightedAdjacencyMatrix` + `EdgeWeight[g]` only, no weighted
shortest-path. See Unresolved above.
