---
source: src/assoc.c
---
**Algorithm.** `builtin_keysort` with one argument copies the entries and `qsort`s
them with `rule_key_cmp`, which orders by `expr_compare` of the keys. Keys are
distinct, so this is a total order and the sort is well defined. The two-argument
`KeySort[assoc, p]` orders by a user ordering function `p`: it extracts the keys and
calls `Ordering[keys, All, p]`, then reindexes the entries by that permutation —
sharing `Sort`'s merge sort, so ties and a non-boolean `p` behave exactly as `Sort`
would.

**Data structures.** An `Expr**` array of entry copies; the system `qsort` (1-arg) or
the `Ordering` permutation vector (2-arg) rebuilt into a fresh `Association`.

**Complexity / limits.** O(n log n). When `p` is undecidable on the keys (e.g.
`Greater` on bare symbols) `Ordering` leaves them in place, so the result keeps the
original order rather than erroring.
