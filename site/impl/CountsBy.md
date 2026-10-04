---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing)."
source: src/assoc.c
---
**Algorithm.** `builtin_countsby` returns `<|f[x] -> count|>`, tallying the
elements by the value of `f` applied to each. It evaluates `f[x]` once per
element (via the `apply1` helper) and keys the result into a `KeyIndex`: a new
`f`-value is recorded with count 1, a repeat increments the stored counter. The
distinct `f`-values and counts become the association's rules in first-appearance
order. An `Association` argument is tallied over its values by `f`
(`assoc_apply_over_values`).

**Data structures.** A `KeyIndex` open-addressing hash set over the *owned*
`f`-value keys, paired with an `int64_t` count array; the keys are adopted into
the result rules so no second copy is made.

**Complexity / limits.** `O(n)` evaluations of `f` plus `O(n)` hashing. Returns
`NULL` unless the argument is a `List` or association. There is no packed fast
path — `f` is an arbitrary symbolic function applied per element — so a packed
argument is materialised before counting.
