---
source: src/funcprog.c
---
`builtin_nestlist` is `nest_impl(res, true)`: it returns the full iteration
history `{expr, f[expr], ..., Nest[f,expr,n]}`. It seeds an `ExprBuf` with a copy
of `expr` and runs the shared `iter_run` driver with `nest_step` (each step is
`apply_unary(f, last)` = build `f[last]` and `eval_and_free`). With
`as_list=true`, `ebuf_finalize` wraps the entire kept history in a `List` head.
Shares all machinery with `Nest`; `n` must be a non-negative integer.
