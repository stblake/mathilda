---
source: src/part.c
---
`builtin_extract` (in `src/part.c`) reads `Extract[expr, pos]` (optionally with a held wrapper head as a third argument). If `pos` is a list-of-positions (a `List` whose elements are themselves `List`s), it extracts each position via the helper `extract_single` and returns the results in a `List`; otherwise it treats `pos` as a single position path. The 1-arg operator form returns a `Function[Extract[#, pos]]` closure.
