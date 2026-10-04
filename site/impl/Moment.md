---
references:
  - "M. G. Kendall, A. Stuart and J. K. Ord, *The Advanced Theory of Statistics*, Vol. 1 (Griffin, 1987), ch. 3 (moments)."
source: src/stats/moment.c
---
**Algorithm.** `builtin_moment` takes `Moment[data, order]`. An `Association`
argument is handled over its values (`assoc_apply_over_values`); an `NDArray` /
packed-array argument takes the buffer fast path `ndred_moment`. Otherwise `data`
must be a non-empty `List` and the dispatch is by the shape of `order`:

1. **List order `{r1, …, rm}`** (`mom_multivariate`) — the multivariate mixed raw
   moment. For each first-axis slice `sub` it forms `Times @@ sub^rvec` (the
   `Listable` `Power` raises block `j` to `r_j` along the second axis, and
   `Times @@` takes the product over that axis), then `Mean` over the slices. The
   slice's second-axis length must equal `Length[rvec]`.
2. **Scalar order over a flat vector** — fast real path: when every element is a
   real number, at least one is inexact, and `order` is a non-negative machine
   integer, a tight C loop accumulates `Σ xᵢ^r` with integer power-by-squaring
   (`mom_ipow`) and returns `acc / n` as a machine `Real`.
3. **Everything else** (`mom_symbolic`) — `Mean[data^r]`. Because the raw moment
   has no mean to subtract, `Power` being `Listable` makes `data^r` thread
   element-wise at *every* rank — so a single expression gives the vector's scalar
   moment, the matrix's columnwise vector, or the array's columnwise array, and no
   separate columnwise routine is needed (unlike `CentralMoment`). This reuses
   `Mean`'s exact-rational, symbolic, and MPFR paths.

`Moment[data, 1]` is `Mean[data]` and `Moment[data, 0]` is `1`, both falling out
of the same definition.

**Data structures.** `Expr` trees driven through `eval_and_free`; the multivariate
case holds an `Expr**` array of `n` per-slice terms before wrapping them in a
`List` for the outer `Mean`. The fast path is a single `double` accumulator.

**Complexity / limits.** `O(n log r)` for the machine-double loop (power by
squaring per element); the symbolic path is `O(n)` tree builds plus the
evaluator's `Mean` cost. Registered as `packed_aware` in `src/pack.c` but **not**
`INT64_OK` — an integer sample's raw moment is a `Rational` no `int64` slot holds,
so an integer buffer degrades to the exact `List` path, exactly like `Variance`.
Lowers inside `Compile[]` at a real vector with integer order, so it participates
in auto-compilation.
