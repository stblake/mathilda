# ConvexHullRegion

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ConvexHullRegion[{{x1, y1}, ...}] gives the convex hull of a set of 2D points: a Polygon with the hull vertices in counterclockwise order, a Line for collinear input, or a Point for a single point.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= ConvexHullRegion[{{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}, {1, 1}}]
Out[1]= Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}]

In[2]:= ConvexHullRegion[{{0, 0}, {1, 1}, {2, 2}, {3, 3}}]
Out[2]= Line[{{0, 0}, {3, 3}}]

In[3]:= ConvexHullRegion[{{1, 2}}]
Out[3]= Point[{1, 2}]

In[4]:= Area[ConvexHullRegion[{{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}, {1, 1}}]]
Out[4]= 4
```

### Applications (4)

Interior and collinear points are dropped

```mathematica
In[5]:= ConvexHullRegion[{{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}, {1, 1}}]
Out[5]= Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}]
```

Collinear input collapses to a Line

```mathematica
In[6]:= ConvexHullRegion[{{0, 0}, {1, 1}, {2, 2}, {3, 3}}]
Out[6]= Line[{{0, 0}, {3, 3}}]
```

A single point gives a Point

```mathematica
In[7]:= ConvexHullRegion[{{1, 2}}]
Out[7]= Point[{1, 2}]
```

Composes with the other region heads

```mathematica
In[8]:= Area[ConvexHullRegion[{{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}, {1, 1}}]]
Out[8]= 4
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

- A. M. Andrew, *Another efficient algorithm for convex hulls in two dimensions*, Inform. Process. Lett. **9** (1979) 216-219.
- Source: [`src/geometry.c`](https://github.com/stblake/mathilda/blob/main/src/geometry.c)
- Specification: [`docs/spec/builtins/geometry.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/geometry.md)
- Tests: [`tests/test_geometry.c`](https://github.com/stblake/mathilda/blob/main/tests/test_geometry.c)

## Notes & additional examples

### Notes

`ConvexHullRegion[{{x1, y1}, …}]` is the convex hull of a 2D point set, computed by Andrew's
monotone chain. Duplicate, interior, and collinear-middle points are removed, and the vertices
come out counterclockwise starting from the lexicographic minimum — the Wolfram Language's own
vertex order.

The result head follows the hull's dimension: a `Polygon` for a genuine 2D hull, a `Line` for
collinear input, and a `Point` for a single distinct point. It returns a bare
`Polygon[{verts}]` — unlike Wolfram 12+, which attaches a cell-spec second argument. Integer/
rational input is processed in exact arithmetic (so the sign tests are exact); a `Real`
coordinate switches to the machine path. A point coordinate too large for a `double` makes the
machine path decline rather than emit a degenerate region. `Protected`, not `Listable`.
