---
source: src/geometry.c
---
**Algorithm.** `builtin_perimeter` unwraps `Polygon[pts]` with `geom_polygon_points`, reads
and dedups the vertices (`geom_read_points` + `geom_dedup_consecutive`, consecutive and
wrap-around repeats collapsed), and returns `Undefined` for fewer than 3 distinct vertices.
The perimeter is the sum of edge lengths including the closing edge,
`Sum_i Sqrt[(x_{i+1}-x_i)^2 + (y_{i+1}-y_i)^2]`:

- **exact path** — each squared edge length is computed in GMP (`mpq_sub`/`mpq_mul`/
  `mpq_add`), wrapped as a `Sqrt[...]` `Expr`, and the `Plus` of those terms is handed to the
  evaluator, so radical canonicalisation does the rest: `Perimeter[Polygon[{{0,0},{1,0},
  {0,1}}]]` becomes `2 + Sqrt[2]` and `Perimeter[Polygon[{{0,0},{1/2,0},{0,1}}]]` becomes
  `3/2 + 1/2 Sqrt[5]`.
- **machine path** — a running `double` sum of `hypot(dx, dy)`, returned as a `Real`.

The exact/machine choice is made once by `geom_read_points` from the coordinate types (WL
contagion: any `Real` coordinate makes the whole result machine).

**Data structures.** A `GeomPoints` with parallel `mpq_t qx/qy` and `double dx/dy`. The exact
path allocates one `Sqrt[...]` `Expr` per edge into a `terms` array, builds a single `Plus`,
and lets `eval_and_free` simplify it (an OOM mid-build frees the partial terms and declines).

**Complexity / limits.** `O(n)` edges, plus the evaluator's radical-simplification cost on the
exact path and GMP arithmetic on each squared length. Simple polygons only; `Protected`, not
`Listable`. A coordinate that overflows a `double` on a machine path makes `Perimeter`
decline.
