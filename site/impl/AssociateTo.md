---
source: src/assoc.c
---
**Algorithm.** `builtin_associate_to` is the in-place update for associations,
mirroring `AppendTo`. It is `HoldFirst`, so the first argument is the unevaluated
symbol; the builtin evaluates it to read the symbol's current association,
gathers the existing entries together with the new rule(s) (a single
`key -> val` or a `List` of rules), re-canonicalises the whole set with
`assoc_from_rules`, assigns the updated association back to the symbol via
`symtab_add_own_value`, and returns it.

**Data structures.** A flat `Expr*` array of the base entries plus the new ones
is passed to `assoc_from_rules`, whose transient `KeyIndex` hash set applies the
*first position, last value* de-duplication — so associating an existing key
overwrites its value in place while a new key is appended at the end.

**Complexity / limits.** `O(n + k)` for `n` existing entries and `k` new rules.
Returns `NULL` (unevaluated) if the first argument is not a symbol, if its
current value is not an association, or if any supplied entry is not a 2-arg
rule.
