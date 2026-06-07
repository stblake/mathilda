---
source: src/calculus/series.c
---
**Data structures.** `SeriesData[x, x0, {a0, ..., a_{k-1}}, nmin, nmax, den]` is
the data head representing a truncated power series produced by Series. The i-th
coefficient `a_i` multiplies `(x - x0)^((nmin + i)/den)`, and the `O[x - x0]^(nmax/den)`
term captures the dropped higher-order tail. The integer `den` (>= 1) is the
common denominator of the exponents, so Laurent (`nmin < 0`) and Puiseux
(`den > 1`, fractional exponents) series are both representable. It carries only
`ATTR_PROTECTED` — there is no `builtin_seriesdata` handler; it is an inert
container constructed by `so_to_expr` from the internal `SeriesObj` and consumed
by `Normal` (which drops the O-term and rebuilds the explicit `Plus` of powers)
and by the printer. The same fields mirror the in-memory `SeriesObj` struct
(`x`, `x0`, owned coefficient array, `nmin`, `order`, `den`) used during
computation in `series_expand`.
