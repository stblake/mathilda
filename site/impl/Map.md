---
source: src/funcprog.c
---
**Algorithm.** `builtin_map` applies `f` to subexpressions of `expr` at the
levels selected by an optional level-spec (default `{1,1}`, the immediate
arguments). The recursion `map_at_level` works **bottom-up**: for an
`EXPR_FUNCTION` it first rebuilds the node by mapping into each argument (and the
head too when `Heads -> True`), then — if the node's current level is within
`[spec.min, spec.max]` (negative levels measured against `get_depth`) — wraps it
in `f[...]` and calls `evaluate()`. Atoms are copied, and tested for membership
of the level range only by their depth.

**Level / option parsing.** `parse_level_spec` reads an integer `n`, `{n}`,
`{m,n}`, or `Infinity`; `parse_options` reads a trailing `Heads -> True`. A
`Rule`-headed third argument is treated as an option rather than a level-spec.

**Data structures.** Pure `Expr`-tree traversal; new nodes built with
`expr_new_function`. Map, MapAll, and MapAt all share this module and the
`LevelSpec { min, max, heads }` struct.
