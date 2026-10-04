---
source: src/assoc.c
---
**Algorithm.** `builtin_keyexistsq` returns `True` iff the literal key is present
in the association (or bare list of rules). It is a single `assoc_lookup_value`
probe: a non-`NULL` value means the key exists. Unlike `KeyMemberQ`/`KeyFreeQ`,
the second argument is treated as a *literal* key, not a pattern — so
`KeyExistsQ[a, _]` looks for the key spelled `_`, not "any key".

**Data structures.** `assoc_lookup_value` uses the association's cached
`AssocIndex` when present (built lazily on the first single-key read), giving an
`O(1)` open-addressing hash probe; otherwise it falls back to an `O(n)`
`assoc_scan`. Keys are compared with `expr_eq`.

**Complexity / limits.** `O(1)` amortised once the key index exists, else `O(n)`.
Returns `NULL` (unevaluated) unless the first argument is an association or a list
of rules; it is a `*Q` predicate and otherwise always answers a Boolean.
