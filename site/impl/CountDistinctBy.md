---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing)."
source: src/assoc_ops.c
---
**Algorithm.** `builtin_countdistinctby` returns the number of distinct values of
`f[element]`. It shares the `count_distinct` core with `CountDistinct`, passing
the second argument as the key function: `f` is evaluated once per element
(`eval_call1`) and the result added to an `ExprSet`; the final set size is the
answer.

**Data structures.** An `ExprSet` hash set over the owned `f`-values
(`expr_hash`/`expr_eq`), plus an array that owns each `f[element]` result for the
duration of the pass. A visible `NDArray` argument is de-listed first
(`ops_delist_visible`).

**Complexity / limits.** `O(n)` evaluations of `f` plus `O(n)` hashing. Requires
a non-atomic first argument (`CountDistinctBy::normal` otherwise, via the message
funnel); `f` is arbitrary, so no packed fast path applies.
