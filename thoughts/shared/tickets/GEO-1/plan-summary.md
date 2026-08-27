---
ticket: GEO-1
created: 2026-08-27
type: plan
lifecycle: active
status: draft
full_plan: thoughts/shared/tickets/GEO-1/plan.md
---

# Plan Summary: Core computational-geometry builtins

**Full plan (appendix)**: `thoughts/shared/tickets/GEO-1/plan.md`

## Recommendation
Add `src/geometry.c` with five builtins — `Area`, `Perimeter`, `RegionCentroid`, `RegionMember` (2D `Polygon`), `ConvexHullRegion` (→ `Polygon`/`Line`/`Point`) — exact via GMP `mpq_t`, machine via doubles, semantics locked to live-verified Wolfram outputs (24 acceptance rows). Three phases: module+wiring, tests (ctest binary), docs+verification ladder.

## Options Considered
1. Five-head `ConvexHullRegion` slice (chosen) — no new object types; hull output feeds the measurement heads.
2. `ConvexHullMesh` + mesh regions — needs a whole object subsystem; rejected.
3. Reuse renderer shoelace — `USE_GRAPHICS`-gated, double-only; rejected.

## Decisions
- New module; `Polygon` stays an inert tag; heads pattern-match on it.
- GMP `mpq_t` exact path (BigInt/Rational-safe), machine-double contagion when any Real present.
- Decline (`return NULL`) for anything out of scope — unevaluated beats wrong.
- Attributes `Protected | ReadProtected` (WL-verified); no Listable.
- `RegionCentroid` declines on zero-area polygons (documented deviation).

## Non-goals
Self-intersecting polygons; 3D; other region heads as measurement args; `Polygon` holes; WL's hull cell-spec second arg; packed/Compile fast paths for these structural heads; the `GeometricRegion` framework.

## Open Questions

### Unresolved
_None._

### Resolved
- [x] Zero-area centroid → decline (deviation, documented). (assumed — beta test)
- [x] Hull vertex order must match WL — monotone chain CCW from lexicographic min, verified per-AC. (assumed — beta test)
- [x] Exact/machine mixing → machine contagion, as WL. (assumed — beta test)

## Requires Approval
No human present (beta run); plan stays `status: draft`. Reviewer to confirm: zero-area-centroid deviation, bare-`Polygon` hull output, simple-polygon-only limitation.

## Architecture Impact
- New services introduced: none
- APIs changed: none (five new heads; no existing head's behavior changes)
- Data crossing a service boundary: none
- New external dependency: none (GMP already unconditional)
- Deployment topology change: none

## Subsystems & Dependencies
- Subsystems touched: none declared (`thoughts/shared/subsystems/` does not exist in this repo; lookup would be empty)
- Interdependencies surfaced: geometry → graphics: none at code level; `Polygon`/`Line`/`Point` symbols shared as inert tags only.
