---
source: src/part.c
---
**Algorithm.** `builtin_most` returns a copy of the input with its last element dropped: it
copies args `0 .. n−2` into a new function node with the same head. Returns `NULL` (unevaluated)
for atoms or empty expressions.
