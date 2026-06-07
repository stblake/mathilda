---
source: src/stats.c
---
**Algorithm.** `builtin_rootmeansquare` computes `Sqrt[Mean[x^2]]`. It reduces matrices column-wise via `apply_columnwise` and requires a `List`. For data containing a real, it sums squares in `double` and returns `expr_new_real(sqrt(sum_sq/n))`. For exact/symbolic data it builds `Plus[x_i^2 ...]`, then carefully distributes the square root to keep results exact: if the summed result is non-numeric and `n` is a perfect square it factors out `1/Sqrt[n]`; if the mean-square is a rational with a perfect-square denominator it pulls that root out before applying `Power[..., 1/2]`. Otherwise it returns `Power[meanSq, 1/2]` for the evaluator to simplify. `ATTR_PROTECTED`.
