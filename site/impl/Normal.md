---
source: src/calculus/series.c
---
**Algorithm.** `builtin_normal` converts a `SeriesData[x, x0, {a0,...,a_{k-1}},
nmin, nmax, den]` into an ordinary polynomial by dropping the O-term. It builds
the base `(x - x0)` (`series_build_xmx0`), then for each non-zero coefficient `a_i`
forms the term `a_i (x - x0)^((nmin+i)/den)` — using an integer exponent when
`den == 1`, otherwise `Rational[num, den]`, and emitting the coefficient bare when
the exponent is 0 — and sums the terms (`Plus`, or the single term / literal `0`
for degenerate cases), evaluating the result. Any argument that is not a 6-element
`SeriesData` is passed through unchanged (`expr_copy`).

**Data structures.** A direct read of the `SeriesData` arg slots (coefficient
`List`, `nmin`, `den`); no `SeriesObj` is reconstructed.
