---
source: src/list/distance.c
---
**Algorithm.** `builtin_edit_distance` returns the Levenshtein distance — the
fewest single-element insertions, deletions and substitutions turning one
sequence into the other. `dist_levenshtein` runs the standard dynamic program
over two rolling rows rather than the full matrix, so it uses O(min(m, n)) memory
and O(m·n) time. Elements are compared with `expr_eq`, so the one routine serves
both strings (compared character by character) and lists of arbitrary
expressions — `EditDistance[{1, 2, 3}, {1, 3}]` is `1`, as in Mathematica.

**Shape and encoding.** `dist_seq_pair` requires both arguments to be strings or
both to be lists; a mixed pair declines. `dist_seq` either borrows a list's
elements or explodes a string into one `Expr` per *byte* — so a multi-byte UTF-8
character counts as several elements, which matches the ASCII and DNA use cases
and is stated rather than silently assumed.

**Complexity / limits.** O(m·n) time, O(min(m, n)) extra memory; returns `-1`
internally on allocation failure, surfaced as an unevaluated call.
`ATTR_PROTECTED`. See `HammingDistance` for the equal-length positional
analogue.
