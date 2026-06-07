---
source: src/funcprog.c
---
**Algorithm.** `builtin_nestwhile` (`nestwhile_impl(res, false)`) iterates `f`
starting from `expr` while a predicate `test` holds, returning the last value. It
accepts up to six arguments `(f, expr, test, m, max, n)`. `m` controls how many
of the most-recent history entries are passed to `test` (an integer, `All`, or
`{mmin, mmax|Infinity}`; default `1`); `max` caps the number of `f`-applications
(integer or `Infinity`); `n` is post-processing — positive `n` applies `f` that
many extra times, negative `n` drops `|n|` iterates from the end.

The shared `iter_run` driver runs `nestwhile_step`: once at least `m_min` history
entries exist it builds `test[recent...]` from the last `min(count, m_max)`
entries, evaluates it, and halts (without applying `f` again) if the result is
not `True`; otherwise it appends `apply_unary(f, last)`. For an unbounded `max`,
`iter_run` is given `ITER_SAFETY_CAP` (1,000,000) as a safety limit, returning
`NULL` if exceeded. Malformed argument specs return `NULL` (unevaluated).
Post-processing reuses `nest_step` (positive `n`) or `ebuf_truncate` (negative
`n`).

**Data structures.** Growable `ExprBuf` history of owned `Expr*`; the trailing
window passed to `test` is a freshly copied argument array per step.
