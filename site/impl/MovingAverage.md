---
source: src/stats.c
---
**Algorithm.** `builtin_moving_average` takes `(list, r)` where `r` is a positive integer window (`EXPR_INTEGER`/`EXPR_BIGINT`) or a `List` of weights. Output length is `n - r + 1`, and the call stays unevaluated unless `1 <= r <= n`. The unweighted form slides a window of `r` elements, builds a sublist, and delegates to `Mean` per window — so it inherits Mean's exact-rational / real / symbolic handling. The weighted form computes `wsum = Plus[w_k]`, the coefficients `w_k / wsum`, and for each window emits `Plus[Times[coef_k, x_{i+k}], ...]`, letting the evaluator simplify. All intermediate trees are built and reduced with `eval_and_free`. `ATTR_PROTECTED`.
