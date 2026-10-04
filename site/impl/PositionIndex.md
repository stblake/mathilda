---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), ch. 11 (open-addressing hash tables)."
source: src/assoc.c
---
**Algorithm.** `builtin_positionindex` builds `<|value -> {positions}|>` in one hash
pass over the list. Each element is probed in a `KeyIndex`; a new element registers a
growable bucket, and every occurrence appends its 1-based position
(`expr_new_integer(i + 1)`). Distinct values keep first-appearance order. The
`assoc_ops_init` wrapper adds `PositionIndex[assoc]`, which maps each distinct value
to the list of *keys* at which it occurs.

**Data structures.** A single open-addressing `KeyIndex` over the distinct elements,
with a per-element growable `Expr**` bucket of integer positions (doubling on growth).

**Complexity / limits.** O(n), amortised O(1) per element. The result is a
re-canonicalised `Association`; because it indexes by value equality, unhashable or
symbolic elements are grouped by `expr_eq` just like concrete ones.
