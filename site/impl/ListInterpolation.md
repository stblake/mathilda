---
references:
  - "C. de Boor, *A Practical Guide to Splines*, rev. ed. (Springer, 2001) — piecewise-polynomial interpolation on a grid."
source: src/interp.c
---
**Algorithm.** `builtin_listinterpolation` is the value-only companion to
`Interpolation`: it interpolates a rectangular array of *values* laid out on a
regular grid and returns an `InterpolatingFunction`. It reads the array's
rectangular shape along a representative spine (`listinterp_shape`, nesting depth =
dimensionality), synthesises the abscissae for each axis (`listinterp_axis` — plain
integer positions `1..n` by default, `n` equally-spaced points from a `{xmin,
xmax}` interval, or an explicit list of positions), stitches the coordinates and
values into the `{{coord, val}, ...}` table that `Interpolation` consumes, and
hands it to the shared `builtin_interpolation_impl`. All numerics, MPFR handling,
options, and the vectorised `InterpolatingFunction[...]` object are therefore
shared verbatim — this is purely a front-end turning "values on a grid" into
"value at abscissa".

Options pass straight through: `InterpolationOrder -> n` sets the
piecewise-polynomial degree (default 3; 0 constant, 1 linear), `Method ->
"Spline" | "Hermite"` picks the scheme, and `PeriodicInterpolation -> True` builds
a periodic interpolant. `List`, packed, and `NDArray` value tensors all feed the
same path, at machine or arbitrary (MPFR) precision matching the data.

**Data structures.** A `shape[]` vector (capped at `LISTINTERP_MAXDIM = 16`
dimensions) and per-axis abscissa `Expr` arrays — exact endpoint `Expr`s copied so
the object's domain prints exactly, interior nodes machine reals — feeding the
`InterpolatingFunction` representation the shared engine builds.

**Complexity / limits.** Dominated by the underlying `Interpolation` build (per
axis, piecewise fits across the grid). The array must be rectangular (enforced
node-by-node by `listinterp_emit`), non-empty on every axis, and at most 16
dimensions; a domain spec must give one `{min, max}` pair (or position list) per
dimension.
