---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing)."
source: src/assoc.c
---
**Algorithm.** `builtin_associationthread` accepts either the two-argument form
`AssociationThread[keys, values]` or the single rule form
`AssociationThread[keys -> values]`, both requiring two equal-length `List`s. It
zips them into `Rule[key_i, value_i]` nodes and canonicalises the result through
`assoc_from_rules`, so a repeated key keeps its first position and takes its last
value (the same de-duplication rule as `Association`).

**Data structures.** The temporary rule array is handed to `assoc_from_rules`,
which drives a transient `KeyIndex` open-addressing hash set over the keys to
de-duplicate in one pass; the rules are deep-copied into the result and the
temporaries freed.

**Complexity / limits.** `O(n)` amortised. Returns `NULL` (unevaluated) unless
both arguments are `List`s of the same length, matching Wolfram's shape
requirement.
