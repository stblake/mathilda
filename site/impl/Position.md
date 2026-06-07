---
source: src/patterns.c
---
**Algorithm.** `builtin_position` (`src/patterns.c`) returns the index paths of every subexpression matching the pattern. It defaults to levels `0..Infinity` and `Heads -> True` (so heads are searched and positions can contain `0`); the level-spec, `Heads` option, and optional result limit `n` are parsed as in `Cases`. The worker `do_position_at_level` carries a running `int64_t* current_path` (the accumulated index trail) and recurses depth-first pre-order: descending into the head appends index `0`, descending into argument `i` appends `i+1`. At each in-range node it calls `match(e, pattern, env)`; on success it materialises the current path as a `List` of integers and appends it to the results. Negative level-specs are resolved against `get_expr_depth_patterns`; collection stops at `max_results`. `Position[pat]` with one argument returns the operator form `Function[Position[#1, pat]]`.

**Data structures.** A reused `int64_t` path array (reallocated one deeper per recursion level), each match snapshotted into an `Expr` `List` of `EXPR_INTEGER`s; a growable `Expr**` results buffer wrapped into the final `List`; one `MatchEnv` per node.
