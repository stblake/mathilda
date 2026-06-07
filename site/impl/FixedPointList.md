---
source: src/funcprog.c
---
**Algorithm.** `builtin_fixedpointlist` (`fixedpoint_impl(res, true)`) iterates
`f` from `expr` until two successive results coincide, returning the entire
history `{expr, f[expr], ..., fixedPoint, fixedPoint}`. Comparison defaults to
`SameQ` (`expr_eq`); a `SameTest -> s` option swaps in `apply_binary(s, last,
next)`. An optional integer `n` bounds the number of applications; `parse_fp_opts`
parses both, rejecting duplicates or malformed specs (returns `NULL`).

It seeds an `ExprBuf` and runs the shared `iter_run` driver with
`fixedpoint_step`: each step computes `next = apply_unary(f, last)` and, if `next`
equals `last` under the chosen same-test, returns `ITER_STEP_HALT_ADD` (push the
final value and stop) — so the fixed point appears as the last two list elements.
Unbounded runs use `ITER_SAFETY_CAP` (1,000,000). Unlike the scalar `FixedPoint`,
the *List* form does **not** intercept `Throw`/`Abort`/`Quit`/`Return` from `f`
(`propagate_throw=false`); such heads are recorded as ordinary history values.

**Data structures.** Growable `ExprBuf` of owned `Expr*`; `ebuf_finalize` wraps
the kept history under a `List` head.
