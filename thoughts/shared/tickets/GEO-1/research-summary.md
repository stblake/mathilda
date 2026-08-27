---
ticket: GEO-1
created: 2026-08-27T17:41:03-0400
researcher: msollami
topic: "Add core computational-geometry builtins (Area, Perimeter, RegionCentroid, RegionMember on Polygon; ConvexHullRegion)"
type: research
lifecycle: active
full_research: thoughts/shared/tickets/GEO-1/research.md
---

# Research Summary: Core computational-geometry builtins

**Full research (appendix)**: `thoughts/shared/tickets/GEO-1/research.md`

## Recommendation
Add a new `src/geometry.c` module registering five builtins — `Area`, `Perimeter`, `RegionCentroid`, `RegionMember` (2D simple `Polygon`) and `ConvexHullRegion` (returns `Polygon`) — with an exact GMP-`mpq_t` path for exact coordinates and a machine-double path otherwise, semantics pinned to 12 live-verified Wolfram kernel outputs. The codebase is greenfield for geometry (Polygon is an inert rendering tag) but every required mechanism exists and has clean exemplars.

## Options Considered
1. Five-builtin Polygon slice (chosen) — smallest coherent WL-12-style vertical slice; no new object types needed.
2. `ConvexHullMesh` + MeshRegion objects — faithful to WL 12.0 naming but requires a whole mesh-object subsystem; far beyond one ticket.
3. Reuse renderer's `polygon_signed_area` — rejected: double-only, private to `USE_GRAPHICS`-gated code that CI builds without.

## Decision Criteria
- Coherence: hull produces a Polygon that the four region measures consume — one closed loop of functionality.
- Exactness: WL returns exact `1/4`, `2 + Sqrt[2]`, `{1/3, 1/3}` for exact input; mathilda's `make_rational_mpz` + symbolic `Sqrt`/`Plus` evaluation make this achievable, and probes confirm it.
- CI safety: no dependence on optional libs (graphics/FLINT/MPFR); GMP is unconditional.

## Open Questions

### Unresolved
_None._

### Resolved
- [x] Slice, exact/machine duality, simple-polygon-only, 2D-only, packed/NDArray posture — all recorded with rationale in the full research doc (each "(assumed — beta test)").

## Requires Approval
No human present (beta run). Reviewer should confirm: five-builtin slice; simple-polygon limitation; `ConvexHullRegion` returning bare `Polygon[{verts}]` (WL adds a cell-spec second arg); `Area` of degenerate polygon → `Undefined`.
