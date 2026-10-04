---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing)."
source: src/assoc_ops.c
---
**Algorithm.** `builtin_countdistinct` returns the number of distinct elements of
a non-atomic expression (the distinct *values*, for an association). It is the
shared `count_distinct` core with a `NULL` function: every element is added to an
`ExprSet` and the final set size is returned as an integer. One hash pass — no
sorting, no `DeleteDuplicates` round-trip.

**Data structures.** `ExprSet`, a hash set over `Expr*` keyed by `expr_hash` and
compared with `expr_eq`, sized once to the element count. A visible `NDArray`
argument is de-listed up front (`ops_delist_visible`) so the count runs over
ordinary elements.

**Complexity / limits.** `O(n)` amortised. A non-atomic first argument is
required; otherwise `CountDistinct::normal` is emitted (through the message
funnel) and the call is left unevaluated.
