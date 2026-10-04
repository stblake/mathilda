---
source: src/assoc.c
---
**Algorithm.** `builtin_keytake` is `key_drop_take(res, take=true)`. For a single
association it calls `assoc_key_select(assoc, karg, take=true)`, which keeps only the
listed keys while preserving the association's order; `KeyDrop` is the same routine
with `take=false`. The key argument may be one key or a `List` of keys. When the
first argument is a non-empty list of associations the call threads, delegating to
itself per element (the column-of-records form).

**Data structures.** The requested keys are loaded into a transient open-addressing
`KeyIndex` once, so each entry is classified by a single O(1) probe. The result is a
fresh `Association`; the compiled evaluator (B3) calls `assoc_key_select` directly,
with no call-node round-trip.

**Complexity / limits.** O(n + k) for `n` entries and `k` requested keys. A first
argument that is neither an association nor a threadable list leaves the call
unevaluated.
