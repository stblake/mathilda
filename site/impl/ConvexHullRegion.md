---
references:
  - "A. M. Andrew, *Another efficient algorithm for convex hulls in two dimensions*, Inform. Process. Lett. **9** (1979) 216-219."
source: src/geometry.c
---
**Algorithm.** `builtin_convex_hull_region` reads a 2D point set with `geom_read_points` and
computes the convex hull by **Andrew's monotone chain**. The point indices are sorted
lexicographically (`x`, then `y`) with `qsort` — `geom_idx_cmp_q` on the exact path,
`geom_idx_cmp_d` on the machine path (there is no `qsort_r` in C99, so the comparison reads
the point set through a file-scope `g_sort_pts` pointer set for the call). Duplicates become
adjacent after the sort and are removed, then the lower and upper chains are built: a vertex
is popped while the last turn is a right turn or collinear (`geom_turn` cross-sign `<= 0`),
which drops interior and collinear-middle points. The hull comes out counterclockwise starting
from the lexicographic minimum — WL's own vertex order. The result head follows the hull size:
`Point[{x, y}]` for one distinct point, `Line[{p1, p2}]` for two (collinear input), and
`Polygon[{verts}]` otherwise.

Both chain loops run the cross-product predicate `O(n)` times, so the exact path uses one
`GeomCrossTmp` block of six `mpq_t` scratch values initialised once per call (per-call
`mpq_init`/`mpq_clear` inside the predicate dominated the loop before this).

**Data structures.** A `GeomPoints` with parallel `mpq_t`/`double` coordinate arrays; an
`idx` array of point indices (sorted and deduped in place) and a `hull` stack of size
`2n + 1`. Hull vertices are emitted as two-element `List`s (`geom_list2`), exact coordinates
via `make_rational_mpz` and machine ones as `Real`s.

**Complexity / limits.** `O(n log n)` for the sort plus `O(n)` for the chain construction;
exact arithmetic on integer/rational input, `double` otherwise. `Protected`, not `Listable`.
Returns a bare `Polygon[{verts}]` without WL 12's cell-spec second argument. A point coordinate
too large for a `double` on the machine path makes it decline rather than emit a degenerate
`Line` with a non-re-parseable infinity.
