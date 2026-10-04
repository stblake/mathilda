---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing)."
source: src/assoc.c
---
**Algorithm.** `builtin_keydrop` is the drop half of `key_drop_take`, the shared
core with `KeyTake`. Given a key or a `List` of keys it builds a `KeyIndex` of
that drop set once, then rebuilds the association keeping exactly the entries
whose key is *absent* from the set — order preserved. A `List` of associations is
threaded element-wise (the column-of-records form), each element delegated back
to the same builtin.

**Data structures.** A transient `KeyIndex` open-addressing hash set over the
requested keys; the surviving entries are deep-copied into a fresh canonical
`Association`. The compiled evaluator calls the lower-level `assoc_key_select`
directly, with no call-node round-trip.

**Complexity / limits.** `O(n + m)` for `n` entries and `m` requested keys,
versus the `O(n·m)` of repeated membership scans. Returns `NULL` unless the first
argument is an association (or a non-empty list of them).
