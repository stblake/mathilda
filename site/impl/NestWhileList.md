---
source: src/funcprog.c
---
`builtin_nestwhilelist` is `nestwhile_impl(res, true)`: same iteration as
`NestWhile` but it returns the full history list `{expr, f[expr], ...}`. It
accepts `(f, expr, test, m, max, n)` with identical semantics — `m` selects how
many recent entries are fed to `test` (integer / `All` / `{mmin,mmax}`, default
`1`), `max` bounds applications, and post-processing `n` appends extra `f`
applications (positive) or trims trailing iterates (negative). The shared
`iter_run` driver runs `nestwhile_step`, which evaluates `test[recent...]` and
halts before the next `apply_unary(f, last)` when the test fails. With
`as_list=true`, `ebuf_finalize` wraps the kept history in a `List`. Unbounded runs
are guarded by `ITER_SAFETY_CAP`; malformed specs return `NULL`.
