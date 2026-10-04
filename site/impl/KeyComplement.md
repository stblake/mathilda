---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing)."
source: src/assoc_ops.c
---
**Algorithm.** `builtin_keycomplement` takes a single `List` of associations
`{a1, a2, ...}` and returns the entries of the first association `a1` whose key
appears in *none* of the later associations — the set-difference of key sets,
carrying `a1`'s values and preserving `a1`'s order. Each candidate key of `a1` is
tested against the other associations with `assoc_lookup_value`, and the entry is
kept only when every other lookup misses.

**Data structures.** The other associations are probed through their cached
`AssocIndex` (open-addressing hash), so each membership test is `O(1)` amortised;
the surviving entries are deep-copied into a fresh `Association`.

**Complexity / limits.** `O(|a1| · m)` lookups for `m` associations, each lookup
`O(1)` amortised. An empty list argument raises `KeyComplement::empt`; a
non-list-of-associations raises `KeyComplement::invar`. Both diagnostics route
through the message funnel (`ops_msg`), so they honour `Quiet[]` and `Check[]`.
