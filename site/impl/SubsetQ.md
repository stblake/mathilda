---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), ch. 11 (hash sets)."
source: src/assoc_ops.c
---
**Algorithm.** `builtin_subsetq` tests whether every element of `b` occurs in `a`,
ignoring multiplicity. It loads all of `a`'s elements into an `ExprSet` (a hash set),
then probes each element of `b`, short-circuiting on the first miss. Lists and
associations mix freely and are compared by element (by *value* for an association);
any other pair of expressions must share a head, or the `SubsetQ::heads` message
fires and the call is left unevaluated.

**Data structures.** An `ExprSet` hash set over `a`'s elements (membership by
`expr_eq`/`expr_hash`); `elem_at` abstracts list vs. association element access.

**Complexity / limits.** O(|a| + |b|) — one pass to build the set, one to probe.
Multiplicity is ignored, so `SubsetQ[{1, 1}, {1}]` and `SubsetQ[{1}, {1, 1}]` are both
`True`.
