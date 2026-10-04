---
source: src/assoc.c
---
**Algorithm.** `builtin_keyfreeq` is the boolean complement of `KeyMemberQ`: both
call `assoc_some_key_matches`, and `KeyFreeQ` negates its verdict. The second
argument is treated as a *pattern*, as in Mathematica. When it is pattern-free
(`key_contains_pattern` finds no `Blank`/`Pattern`/`Alternatives`/… construct) the
test is a single `assoc_lookup_value` probe — the O(1) index lookup. When it does
carry a pattern, the matcher is run against each key in turn (`match` into a fresh
`MatchEnv`) and the first match wins. Either an association or a bare list of rules
is accepted.

**Data structures.** For the literal-key path, the association's persistent
open-addressing hash index (`KeyIndex`, keys compared by `expr_eq` and hashed by
`expr_hash`); for the pattern path, a throwaway `MatchEnv` per key.

**Complexity / limits.** O(1) amortised for a literal key, O(n) for a pattern (one
structural match per entry). Contrast `KeyExistsQ`, which looks its argument up
literally even when it is a `_`.
