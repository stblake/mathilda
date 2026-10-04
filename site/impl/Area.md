---
source: src/geometry.c
---
**Algorithm.** `builtin_area` unwraps `Polygon[pts]` (exactly one part) via
`geom_polygon_points`, reads the vertices with `geom_read_points`, and collapses consecutive
duplicate vertices — including the wrap-around pair — with `geom_dedup_consecutive`, so an
explicitly closed vertex list and interior repeats are handled. A polygon with fewer than 3
distinct vertices is degenerate and returns `Undefined`. The area is the shoelace sum
`|Sum_i (x_i y_{i+1} - x_{i+1} y_i)| / 2`:

- **exact path** (every coordinate `Integer`/`BigInt`/`Rational`) — `geom_signed_area2_q`
  accumulates twice the signed area in a GMP `mpq_t`, which is then made positive, halved
  (`mpq_div_2exp`), and canonicalised to an `Integer` or `Rational` by `make_rational_mpz`. So
  `Area[Polygon[{{0,0},{1,0},{1/2,1/2}}]]` is exactly `1/4`.
- **machine path** (any `Real`/`MPFR` coordinate, WL's contagion) — `geom_signed_area2_d`
  sums in `double` and the result is `fabs(area2)/2` as a `Real`.

`geom_read_points` classifies the whole vertex list in one pass (exact until a `Real` forces
machine; decline on a symbolic or complex coordinate), so the two paths never mix.

**Data structures.** A `GeomPoints` holding parallel `mpq_t qx/qy` (exact) and `double
dx/dy` (always mirrored), with an `exact` flag and a `d_finite` flag. A visible `NDArray` of
points is read directly: an `int64`-typed array keeps the exact path (so it agrees with the
same points as a nested `List`), any other dtype takes the machine path.

**Complexity / limits.** `O(n)` in the number of vertices, plus GMP cost on the exact path.
Simple (non-self-intersecting) polygons only — a self-intersecting vertex list follows
shoelace/even-odd semantics, not WL's enclosed-region model, and is not detected at runtime.
`Area` is `Protected` and deliberately **not** `Listable` (it is structural, not
element-wise); a coordinate too large for a `double` on a machine path makes it decline rather
than compute on an infinity.
