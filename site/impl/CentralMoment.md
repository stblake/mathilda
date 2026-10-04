---
references:
  - "M. G. Kendall, A. Stuart and J. K. Ord, *The Advanced Theory of Statistics*, Vol. 1 (Griffin, 1987), ch. 3 (moments about the mean)."
source: src/stats/central_moment.c
---
**Algorithm.** `builtin_central_moment` takes `CentralMoment[data, order]`, the
`order`-th moment about the mean, `μ̃_r = (1/n) Σ (xᵢ − Mean[data])^r`. An
`Association` is handled over its values; an `NDArray` / packed argument takes the
buffer fast path `ndred_central_moment`. Otherwise `data` must be a non-empty
`List` and dispatch is by shape:

1. **List order `{r1, …, rm}`** (`cm_multivariate`) — the multivariate mixed
   central moment. With `mu = Mean[data]` (rank `k−1`), each first-axis slice
   contributes `Times @@ (sub − mu)^rvec`, then `Mean` over the slices. The slice's
   second-axis length must equal `Length[rvec]`.
2. **Scalar order over an array of depth ≥ 2** (`cm_columnwise`) — reduces the
   first axis columnwise: `mu = Mean[data]`, and each slice gives `(sub − mu)^r`,
   `Mean`-averaged over the slices. The subtraction is done *per slice* because the
   obvious `data − Mean[data]` would thread row-wise, not column-wise.
3. **Scalar order over a flat vector** — fast real path when every element is a
   real number, at least one is inexact, and `order` is a non-negative machine
   integer: one pass computes the mean, a second accumulates `Σ (xᵢ − mean)^r` by
   integer power-by-squaring (`cm_ipow`), returning `acc / n` as a machine `Real`.
4. **Everything else** (`cm_vector_symbolic`) — `Mean[(data − Mean[data])^r]`,
   where the `Listable` `Plus` and `Power` thread over the vector and the outer
   `Mean` divides by `n`, reusing `Mean`'s exact and symbolic paths.

The design mirrors `Variance` with the `n/(n−1)` bias correction removed: it
divides by `n` (not `n−1`), raises to `r` (not a square), needs only `n ≥ 1`, and
applies no `Conjugate` — a central moment is `(x − μ)^r`, not `|x − μ|²`.
`CentralMoment[data, 1]` is `0`.

**Data structures.** `Expr` trees via `eval_and_free`; the columnwise and
multivariate cases hold an `Expr**` of `n` per-slice terms. The fast path uses two
`double` scalar passes.

**Complexity / limits.** `O(n log r)` machine-double, two passes; the symbolic path
is `O(n)` tree builds plus `Mean`. `packed_aware` in `src/pack.c` but not
`INT64_OK` (an integer sample's central moment is a `Rational`, so an integer
buffer degrades to the exact `List` path). Lowers inside `Compile[]` at a real
vector with integer order; participates in auto-compilation.
