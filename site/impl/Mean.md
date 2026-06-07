---
source: src/stats.c
---
**Algorithm.** `builtin_mean` first probes its argument with `MatrixQ`; if true it computes column-wise means via `apply_columnwise` (which is `Map[Mean, Transpose[matrix]]`). Otherwise it requires a `List` (`ListQ`). For a vector of length `n` it dispatches on element kinds: if any element is `EXPR_REAL`, it sums to a `double` and returns `expr_new_real(sum/n)`; if all elements are exact integers/rationals it accumulates the sum in `int64_t` numerator/denominator pairs (reducing by `gcd` each step) and returns `make_rational(sum_n, sum_d * n)`. Anything symbolic falls back to `(1/n) * (Plus @@ data)` built as `Times`/`Apply` nodes and re-evaluated.

**Limits.** The exact-rational accumulator uses fixed `int64_t` arithmetic, so it can overflow for large/many rationals (no GMP promotion in this path). Empty list returns `NULL`.
