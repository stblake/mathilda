---
references:
  - "M. G. Kendall, A. Stuart and J. K. Ord, *The Advanced Theory of Statistics*, Vol. 1 (Griffin, 1987), §3.31 (the standardized third moment)."
source: src/stats/skewness.c
---
**Algorithm.** `Skewness[data]` is the coefficient of skewness, a measure of
asymmetry, equal to `CentralMoment[data, 3] / CentralMoment[data, 2]^(3/2)`.
`builtin_skewness` is a one-line call into the shared body
`stats_standardized_moment(res, 3)` in `src/stats/stats_common.c`, which `Kurtosis`
reuses with `p = 4`. That body:

1. handles an `Association` over its values (`assoc_apply_over_values`);
2. routes an `NDArray` / packed argument to the buffer kernel `ndred_skewness`;
3. otherwise forms `mp = CentralMoment[data, 3]` and `m2 = CentralMoment[data, 2]`
   and, if *either* comes back still wearing a `CentralMoment` head (the data did
   not reduce — a bare symbol, an empty list), returns `NULL` so the caller's head
   stays unevaluated;
4. returns `mp / m2^(3/2)`.

Because `Power` and `Divide` thread, a matrix's columnwise vector of central
moments yields a columnwise vector of skewnesses. Exact input gives exact output,
which is a radical in general (e.g. a `Sqrt[2]` multiple).

**Data structures.** `Expr` trees built and evaluated through `eval_and_free`; the
two central moments are the only intermediates. The buffer path works on the
packed `double` sample inside `ndred_skewness`.

**Complexity / limits.** Two `CentralMoment` reductions, each `O(n)`, plus a
constant-cost combine. `packed_aware` in `src/pack.c` but not `INT64_OK` — an
integer sample's skewness is a radical no machine slot holds, so an integer buffer
degrades to the exact `List` path. Lowers inside `Compile[]` (a `ND_REDS`
delegate) and participates in auto-compilation. A symmetric sample gives exactly
`0`.
