---
source: src/assoc.c
---
**Algorithm.** `builtin_keymemberq` asks whether any key of the association matches
its second argument, which is read as a *pattern*. `assoc_some_key_matches` first
calls `key_contains_pattern`: a pattern-free query is answered by a single
`assoc_lookup_value` probe (O(1) index lookup); a query carrying a pattern
(`Blank`, `Pattern`, `Alternatives`, `PatternTest`, …) is matched structurally
against each key, short-circuiting on the first hit. An association or a bare list
of rules is accepted.

**Data structures.** The association's persistent open-addressing hash index
(`KeyIndex`) for the literal path; a transient `MatchEnv` per key for the pattern
path.

**Complexity / limits.** O(1) amortised for a literal key, O(n) for a pattern.
Unlike `KeyExistsQ` (literal argument), `KeyMemberQ[a, _]` is `True` for any
non-empty `a`.
