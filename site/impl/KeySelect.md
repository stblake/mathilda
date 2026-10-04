---
source: src/assoc.c
---
**Algorithm.** `builtin_keyselect` keeps the entries whose *key* passes a predicate.
For each entry it forms `pred[key]` (`apply1`), evaluates it, and keeps the whole
entry only when the verdict is literally the symbol `True`; any other result drops
the entry. The surviving entries, in their original order, form a fresh
`Association`.

**Data structures.** A single `Expr**` output array sized to the entry count; the
predicate application is a throwaway evaluated per key.

**Complexity / limits.** O(n), one predicate evaluation per entry. The companion
`KeySelect` tests keys; `Select` (and `Discard`) over an association test the values
instead. A non-`True` verdict — including an unevaluated predicate call — counts as
rejection.
