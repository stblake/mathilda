---
source: src/funcprog.c
---
**Algorithm.** `builtin_apply` replaces the head of `expr` with `f` at the
levels given by an optional level-spec (default level `0`, i.e. just the top
head). The work is done by the recursive `apply_at_level`, which descends each
`EXPR_FUNCTION`. A node is "at" the active range when its current depth lies in
`[spec.min, spec.max]`, with negative levels measured against `get_depth(expr)`
(so `-1` etc. count from the leaves). When a node is in range its arguments are
first transformed recursively, a fresh `f[args...]` is built, and `evaluate()`
is called on it so `f`'s attributes take effect; otherwise the original head is
kept (or transformed too when `Heads -> True`).

**Level / option parsing.** The third argument is interpreted by
`parse_level_spec` (handles an integer `n`, `{n}`, `{m,n}`, `Infinity`, and
treats `Automatic`/missing as level `{0,0}`); any trailing `Heads -> True`
option is read by `parse_options`. A `Rule`-headed third argument is recognised
as an option, not a level-spec.

**Data structures.** Operates directly on the `Expr` tagged union; rebuilds
`EXPR_FUNCTION` nodes with `expr_new_function`. Leaves (non-function atoms) are
returned via `expr_copy` unchanged.
