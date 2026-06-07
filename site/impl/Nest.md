---
source: src/funcprog.c
---
`builtin_nest` (via `nest_impl(res, false)`) applies `f` to `expr` exactly `n`
times and returns the final value; `n` must be a non-negative integer or the
call stays unevaluated. It seeds a growable history buffer `ExprBuf` with a copy
of `expr` and drives the shared generic runner `iter_run` with `nest_step`, which
on each step computes `apply_unary(f, last)` (build `f[last]`, `eval_and_free`).
`ebuf_finalize(..., as_list=false)` returns the last history element and frees the
rest. `Nest` and `NestList` are the same `nest_impl`, differing only in the
`as_list` flag.
