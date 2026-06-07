---
source: src/picostrings.c
---
`builtin_stringjoin` gathers all leaf strings into a growable `const char**` array via the recursive helper `collect_strings`, which descends through any `List` wrappers and borrows (does not copy) each `EXPR_STRING`'s `data.string`; any non-string, non-`List` leaf aborts with `NULL` (unevaluated). It then sums the lengths, `malloc`s one buffer of `total_len + 1`, `memcpy`s each fragment in order, and returns a single `EXPR_STRING`. The zero-argument form yields `""`. Registered with `ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED`, so the evaluator flattens nested `StringJoin` (and the `<>` infix operator) before the builtin runs.
