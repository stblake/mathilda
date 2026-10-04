---
source: src/assoc.c
---
**Algorithm.** `builtin_keymap` applies `f` to each key (`apply1(f, key)`), keeps the
value unchanged, and rebuilds the association with `assoc_from_rules`. The
re-canonicalisation matters: `f` may map two distinct keys to the same new key, and
`assoc_from_rules` then applies the usual association rule — a later entry's value
wins for a collided key.

**Data structures.** A fresh `Rule[newkey, value]` array fed to `assoc_from_rules`,
which de-duplicates keys through its own hash index.

**Complexity / limits.** O(n), one `f` evaluation per key, plus the O(n) canonicalise.
Because keys can collide, the result may be shorter than the input. `KeyValueMap`
instead feeds both key and value to `f` and returns a plain `List`.
