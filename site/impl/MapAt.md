---
source: src/funcprog.c
---
**Algorithm.** `builtin_map_at` applies `f` at explicit positions rather than at
levels. It first disambiguates a single position from a list of positions: the
argument is treated as *multiple* paths only when it is a non-empty `List` whose
first element is itself a `List`. A single path is then a position vector
`{i1, i2, ...}` (or a bare index), and the recursive `mapat_at_path` walks it:
when the path is exhausted it applies `f` to the targeted node
(`mapat_apply_f`, which builds `f[node]` and calls `evaluate()`); otherwise it
rebuilds the current `EXPR_FUNCTION` with the chosen child replaced by the
recursive result. A path step may be a positive/negative integer (negatives
count from the end, `0` targets the head), the symbol `All` (apply to every
child at that level), or a `Span[a, b]` / `Span[a, b, step]` range. Out-of-range
indices are silently ignored, matching Mathematica's permissive behaviour.

For the multiple-positions form the paths are applied **sequentially** to a
running copy of the expression, so repeated positions apply `f` more than once.

**Data structures.** Operates on the `Expr` tree; copies each level's argument
array and overwrites only the targeted slot, then rebuilds with
`expr_new_function`.
