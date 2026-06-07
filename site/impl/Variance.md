---
source: src/stats.c
---
**Algorithm.** `builtin_variance` computes the *sample* variance (divisor `n-1`, requiring `n > 1`). It first reduces matrices column-wise via `apply_columnwise`, then requires a `List`. For real-valued data it uses Welford's online algorithm (running mean `m` and sum-of-squares `s`) and returns `expr_new_real(s/(n-1))`. For exact integer/rational data it does the computation in `int64_t` numerator/denominator pairs: it first accumulates the sum (hence the mean), then accumulates `Sum[(x_i - mean)^2]` as exact rationals (reducing by `gcd`), and returns `make_rational(sq_sum_n, sq_sum_d * (n-1))`. The symbolic fallback evaluates `Mean[data]`, forms `Sum[(x - mu) Conjugate[x - mu]]` (so complex/symbolic data gives the Hermitian variance), and divides by `n-1`.

**Limits.** The exact path uses fixed-width `int64_t`, so large rationals can overflow; `n <= 1` returns `NULL`.
