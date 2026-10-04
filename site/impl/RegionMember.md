---
source: src/geometry.c
---
**Algorithm.** `builtin_region_member` takes `RegionMember[Polygon[pts], {x, y}]`, reads and
dedups the polygon vertices, declines for fewer than 3 distinct vertices, and reads the query
point with `geom_read_query`. It is a boundary-inclusive point-in-polygon test by the
crossing-number (ray-casting) method, with an explicit on-edge check:

- **exact path** (`geom_member_q`, used when both polygon and query are exact) — for each edge
  it takes the sign of the cross product `(v_{i+1}-v_i) × (p-v_i)` with `geom_cross_sign_q`
  (exact GMP). A zero cross with the point inside the edge's bounding box means the point lies
  **on** the boundary → immediately `True`. Otherwise an upward/downward crossing of the
  horizontal ray toggles the inside flag, with the crossing direction matched to the cross
  sign so a vertex is counted once. The vertices themselves count as members.
- **machine path** (`geom_member_d`) — the identical logic in `double`, using **exact IEEE
  comparisons** (a cross product of exactly `0.0` is on-edge; no epsilon tolerance).

It returns the symbol `True` or `False`. An exact polygon whose coordinates do not fit a
`double`, queried with a machine point, **declines** rather than compute NaN cross products
that would read as "collinear" and answer wrongly.

**Data structures.** A `GeomPoints` with parallel `mpq_t`/`double` coordinates and a
`d_finite` flag; a `GeomCrossTmp` block of six `mpq_t` scratch values initialised once and
reused across the whole edge walk (the exact cross-sign predicate sits in an `O(n)` loop, so
per-call `mpq_init`/`mpq_clear` would dominate it).

**Complexity / limits.** `O(n)` edges per query. Simple polygons only; `Protected`, not
`Listable`. Declines on a zero-area/degenerate polygon and on the oversized-coordinate machine
case above. Boundary membership is exact on exact input.
