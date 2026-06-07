---
source: src/funcprog.c
---
`builtin_map_all` is a thin wrapper around the same `map_at_level` traversal used
by `Map`, but with the fixed level-spec `{0, Infinity}` (`min=0`,
`max=1000000`, `heads=false`), i.e. `MapAll[f, expr]` ≡ `Map[f, expr, {0,
Infinity}]`. The bottom-up recursion rebuilds every `EXPR_FUNCTION` from its
mapped children and then wraps each node (including the whole expression at level
0) in `f[...]`, calling `evaluate()` so `f`'s attributes apply. A trailing
`Heads -> True` option is honoured via `parse_options`.
