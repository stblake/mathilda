---
source: src/stats.c
---
**Algorithm.** `builtin_moving_median` takes `(list, r)` with `r` a positive integer window (`EXPR_INTEGER`/`EXPR_BIGINT`); output length is `n - r + 1` and the call is unevaluated unless `1 <= r <= n`. It auto-detects vector vs. matrix mode by whether the first element is a `List`, and validates every leaf with `is_real_numeric` (matrices must be rectangular); invalid data prints `MovingMedian::arg1` and leaves the call unevaluated. It then slides the window, builds an `r`-element sublist, and delegates each window to `Median` (so a matrix window yields a column-wise median vector). Window results are assembled into the output `List`. `ATTR_PROTECTED`.
