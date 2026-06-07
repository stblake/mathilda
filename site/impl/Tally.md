---
source: src/list.c
---
**Algorithm.** `builtin_tally` counts distinct elements, returning `{element, multiplicity}`
pairs in first-occurrence order. With the default sameness test it uses a chained hash table
(`expr_hash` for bucketing, `expr_eq` for equality) for `O(n)` expected counting; with a custom
two-argument test it falls back to an `O(n²)` linear scan, evaluating `test[a, b]` per
comparison. Multiplicities are tracked in a parallel `int64_t` array.
