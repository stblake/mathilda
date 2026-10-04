---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), ch. 11 (open-addressing hash tables)."
source: src/assoc_ops.c
---
**Algorithm.** `builtin_keyintersection` restricts every association in a list to the
keys common to all of them. `as_assoc_array` first normalises the argument (its
elements may be associations, rules or lists of rules). It then scans the first
association's keys in order, keeping a key only when `assoc_lookup_value` finds it in
every other association; those common keys (in the first association's order) become
the key set each input is rebuilt against.

**Data structures.** Each association's persistent hash index backs the O(1)
membership probes; an `ExprBuf` holds the common keys, and one `ExprBuf` per
association accumulates the restricted entries.

**Complexity / limits.** O(k₀ · m) index probes for `m` associations with `k₀` keys in
the first. A malformed argument (not a list of associations or rules) yields
`KeyIntersection::invar` and leaves the call unevaluated. `KeyComplement` is the
sibling that keeps the first association's keys found in *none* of the others.
