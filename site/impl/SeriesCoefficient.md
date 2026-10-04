---
source: src/calculus/series.c
---
**Algorithm.** `builtin_seriescoefficient` implements `SeriesCoefficient[f, {x, x0, k}]` —
the coefficient of `(x − x0)^k` — by reusing `Series`. It first type-checks the expansion
variable with `series_spec_vars_ok` (the same gate `Series` uses: the variable must evaluate
to a symbol), placed *before* the symbolic-index special cases so that, with e.g. `x = 5` in
scope, a spec `{5, 0, n}` cannot spuriously match. Two closed-form general terms are then
recognised directly for a symbolic index `n` at `x0 = 0` and emitted as a `Piecewise`:
`ProductLog[x]` gives `(−n)^(n−1)/n!` for `n ≥ 1`, and `FresnelC`/`FresnelS[x]` give their
`4m+1` / `4m+3` power-class coefficients. Inexact arguments route through
`internal_rationalize_then_numericalize`. The general path parses the spec with
`parse_series_spec`, declines the leading `{x, x0}` two-argument form, expands with
`Series[f, {x, x0, max(k, 0)}]` (`eval_and_free`), and reads the coefficient out of the
resulting `SeriesData`: with minimal power `nmin` and denominator `den`, power `k` sits at
list index `j = k·den − nmin`, so an in-range `j` copies that coefficient and an out-of-range
`j` returns `0`. If `Series` collapses to something free of `x`, the result is that value for
`k = 0` and `0` otherwise.

**Data structures.** `Expr` trees throughout. The intermediate is the `SeriesData` six-tuple
`{var, x0, {c_0, …}, nmin, nmax, den}`; the index arithmetic (`k·den − nmin`) is done in
`long long`. The two symbolic-index cases build a `Piecewise[{{value, cond}}, 0]` `Expr`.
`SeriesCoefficient` is registered `Protected` and, unlike Mathematica folklore, does **not**
hold its arguments — the variable is type-checked instead (see the `series_init` comment),
so a named spec and a `Sequence @@` splice both work.

**Complexity / limits.** Cost is that of one `Series` expansion to order `k` (a full Taylor
pass through the series engine), followed by `O(1)` extraction; the result is exact (rational
or symbolic). Composite results that are a prefactor times a `SeriesData` (e.g. asymptotic
expansions at `Infinity`, `1/x` expansions), and non-integer indices, are left unevaluated
(`NULL`). The symbolic-index general term is produced only for the three heads above; any
other symbolic `n` declines. Being a structural head returning symbolic coefficients, it
carries no NDArray/packed kernel and no `Compile[]` lowering.
