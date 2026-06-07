---
source: src/sum/sum_gosper.c
---
**Algorithm.** `DifferenceDelta[f, i]` computes the forward difference
`(f /. i -> i+1) - f`, the discrete analogue of `D` and the left inverse of
indefinite `Sum`. `builtin_differencedelta` (src/sum/sum_gosper.c) requires the
second argument to be a symbol, substitutes `i -> i+1` into `f` (`shift_var`,
implemented via `ReplaceAll`), subtracts the original `f`, and returns the
`Expand`-ed result. Returns NULL (unevaluated) unless the variable is a symbol.

**Data structures.** Plain `Expr*` tree manipulation built on the existing
`ReplaceAll`, subtraction, and `Expand` builtins (`shift_var`, `sum_sub`,
`sum_eval`). No closed-form machinery — it is a thin structural operator that
lives alongside Gosper's summation because the two are inverse operations.
