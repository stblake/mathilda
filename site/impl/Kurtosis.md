---
references:
  - "M. G. Kendall, A. Stuart and J. K. Ord, *The Advanced Theory of Statistics*, Vol. 1 (Griffin, 1987), §3.31 (the standardized fourth moment)."
source: src/stats/kurtosis.c
---
**Algorithm.** `Kurtosis[data]` is the coefficient of kurtosis — a measure of how
concentrated the distribution is in the peak and tails versus the flanks — equal
to `CentralMoment[data, 4] / CentralMoment[data, 2]^2`. This is **Pearson**
kurtosis, not the excess form (which subtracts `3`). `builtin_kurtosis` is a
one-line call into the shared body `stats_standardized_moment(res, 4)` in
`src/stats/stats_common.c`, the same routine `Skewness` uses with `p = 3`:

1. an `Association` is handled over its values;
2. an `NDArray` / packed argument routes to the buffer kernel `ndred_kurtosis`;
3. otherwise `mp = CentralMoment[data, 4]` and `m2 = CentralMoment[data, 2]` are
   formed, and if either stays a `CentralMoment` head (data not reducible) the
   result is `NULL` and the caller's head is left unevaluated;
4. the result is `mp / m2^(4/2) = mp / m2²`.

`Power` and `Divide` thread, so a matrix yields the columnwise vector of kurtoses.
Exact input gives exact output (a `Rational` here, since both moments enter with
even powers).

**Data structures.** `Expr` trees through `eval_and_free`; the two central moments
are the only intermediates. The buffer path works on the packed `double` sample
inside `ndred_kurtosis`.

**Complexity / limits.** Two `CentralMoment` reductions (`O(n)` each) plus a
constant combine. `packed_aware` in `src/pack.c` but not `INT64_OK` (the integer
answer is a `Rational`, so an integer buffer degrades to the exact `List` path).
Lowers inside `Compile[]` (a `ND_REDS` delegate) and participates in
auto-compilation.
