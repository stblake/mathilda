---
references:
  - "R. J. Hyndman and Y. Fan, *Sample quantiles in statistical packages*, The American Statistician **50** (1996) 361–365."
source: src/stats/quantile.c
---
**Algorithm.** `builtin_quantile` takes `Quantile[data, q]`, `Quantile[data,
{q1, …}]`, or `Quantile[data, q, {{a,b},{c,d}}]`. A visible `NDArray` / packed
argument is materialised to the exact `List` path with `pack_unpack`
(correctness-first — no buffer kernel yet; the audit baselines carry the declared
reason). `data` must be a non-empty `List`.

1. **Matrix** (first element is itself a `List`) — transpose and recurse
   columnwise, carrying `q` and the parameter matrix, so each column is quantiled
   independently.
2. **Validation** — every element must be real-numeric (`Quantile::rectn`
   otherwise, including an already-evaluated `Complex[re, im]` with nonzero `im`).
3. **Parameters** — `parse_param_matrix` validates a `{{a,b},{c,d}}` spec (head and
   shape checked, all four entries real-numeric); the default is `{{0,0},{1,0}}`,
   Hyndman–Fan **Type 1**, the left-continuous inverse CDF.
4. **Sort once** with `pack_eval_plain` (so a large machine-number sort's packed
   result is read through `.args`).
5. Each `q` goes through `quantile_one` → `stats_quantile_point` (shared with
   `Quartiles`, in `stats_common.c`). An exact-irrational `q` (e.g. `1/Sqrt[2]`) is
   read numerically through `N[]`; `q ∉ [0, 1]` raises `Quantile::q100`.

`stats_quantile_point` computes `h = a + (n+b)q`, edge-clamps (`h ≤ 1` → first
element, `h ≥ n` → last), takes `j = Floor[h]` clamped to `[1, n−1]`, `g = h − j`,
and the weight `w = c + d·g`. The upper index is `j` when `g = 0` (the two
neighbours coincide at integer `h`) else `j+1`. At `w = 0` or `w = 1` the element
is **selected** outright. For `w ∈ [0, 1]` it returns the convex combination
`(1−w)·A[j] + w·A[j+1]`; outside `[0, 1]` (which the parameterization permits,
though no standard type uses it) it returns the difference form `A[j] + w·(A[j+1] −
A[j])`. Each form is used only where it is numerically safe — the convex form
overflows on neighbours of opposite huge magnitude, the difference form on `w`
outside the unit interval. All arithmetic runs in the evaluator (exact in, exact
out); `double`s are read only for the clamp, index, and weight-regime decisions.

**Data structures.** `Expr` trees through `eval_and_free`; one sorted copy of the
data; an `Expr**` of per-`q` results for the list form. The four parameters are
kept as fresh `Expr*` copies freed on every exit path.

**Complexity / limits.** `O(n log n)` to sort, then `O(1)` per quantile point.
`Median` deliberately differs from `Quantile[…, 1/2]` on even-length data (`Median`
averages the two central order statistics; the default `Quantile` type selects one).
No `Compile[]` lowering and no packed kernel — a visible `NDArray` is unpacked
first.
