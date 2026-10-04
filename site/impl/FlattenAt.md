---
source: src/list/flatten_at.c
---
**Algorithm.** `builtin_flatten_at` splices the subexpression at one or more
positions into its parent, removing that subexpression's head — it works on any
head, not just `List` (`FlattenAt[f[g[1, 2], g[3, 4]], 1]` is `f[1, 2, g[3, 4]]`).
Position resolution — an integer (negatives counting from the end), a single deep
path `{i, j, ...}`, or a list of paths `{{i1, ...}, {i2, ...}}` — is delegated to
the shared walker `expr_apply_at_positions` (`src/part.c`), the same one `MapAt`
and `ReplaceAt` use, so the one-path/many-paths distinction is exactly theirs; an
out-of-range position makes the walker return `NULL` and `FlattenAt` stay
unevaluated. The splice itself is deferred to the evaluator: the leaf action
`flatten_at_leaf` replaces the addressed `g[a1, ..., ak]` with
`Sequence[a1, ..., ak]`, and `flatten_sequences` (`src/eval.c`) splices that
`Sequence` into the parent on the next pass. Because each targeted slot holds a
single `Sequence` node during the walk, `arg_count` never changes mid-walk, so
several positions resolve against the *original* expression with no index-shift
bookkeeping (unlike `Insert`/`Delete`, which must sort positions descending).

**Data structures.** No bespoke structure — the walker owns the traversal and the
new tree is built by `expr_copy` of the leaf's arguments into the `Sequence`.
`flatten_at_atomic` guards the `Rational`/`Complex` nodes (stored as
`EXPR_FUNCTION` but atomic), so `FlattenAt[{1/2, x}, 1]` never manufactures
`Sequence[1, 2]`.

**Complexity / limits.** Linear in the addressed subexpressions' sizes. A visible
`NDArray` is atomic, so it is first materialised once to a nested `List`
(`ndarray_to_nested_list`) and flattened there; the result is ragged by
construction, so it is deliberately never repacked. A structural head — not on
`pack.c`'s `AWARE` list and with no `Compile[]` lowering. `Protected`.
