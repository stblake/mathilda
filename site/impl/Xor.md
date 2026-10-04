---
source: src/boolean.c
---
**Algorithm.** `Xor` is `Flat, Orderless, OneIdentity, Protected`, so by the time
`builtin_xor` runs the evaluator has already flattened nested `Xor` calls and
sorted the arguments canonically. The builtin makes one pass over the arguments:
each literal `True` flips a parity bit (and is dropped), each `False` is dropped,
and any argument structurally equal (`expr_eq`) to one already kept cancels
against it (`a` Xor `a` = `False`). What survives is the deduplicated non-literal
core; if an odd number of `True`s were consumed, the result is wrapped in `Not`.
`Xor[]` is `False`, `Xor[e]` collapses to `e` by `OneIdentity`, and when nothing
simplified the builtin returns `NULL` to stay symbolic.

**Data structures.** A single `malloc`'d array of *borrowed* argument pointers
holds the surviving terms; cancelled slots are `NULL`'d and then compacted, and
only the final core is deep-copied into the result. The parity is one `int`.

**Complexity / limits.** The duplicate check is pairwise, so the pass is `O(n²)`
in the argument count with `expr_eq` comparisons — fine for the small boolean
expressions this is used on. Simplification is purely structural: it folds literal
Booleans and exact duplicates but does not reason about implications between
distinct symbolic arguments.
