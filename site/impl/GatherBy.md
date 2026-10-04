---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing)."
source: src/assoc.c
---
**Algorithm.** `builtin_gatherby` delegates to `assoc_gather_core`, the single
grouping engine shared with `Gather`. It evaluates `f[x]` once per element, looks
the key up in a `KeyIndex`, and appends the element to that key's buffer —
exactly like `GroupBy`, but the result is the plain list of groups with the group
keys dropped. Groups appear in first-appearance order and keep input order
within each group. Over an `Association` the entries are gathered by `f[value]`
into sub-associations (keys preserved), returned as a list.

**Data structures.** A `KeyIndex` open-addressing hash set over owned keys, with
per-group doubling `Expr**` buffers; the outer result is a `List` of `List`s (or
of sub-associations).

**Complexity / limits.** `O(n)` evaluations of `f` plus `O(n)` hashing. Returns
`NULL` for a non-2-argument call or a non-list/association first argument. `f` is
arbitrary, so no packed fast path applies. `Gather[list]` is the identity-key
special case (`assoc_gather_core(list, NULL)`), which skips the per-element
function application entirely.
