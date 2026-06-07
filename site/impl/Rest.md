---
source: src/part.c
---
**Algorithm.** `builtin_rest` returns a copy of the input with its first element dropped: it
copies args `1 .. n−1` into a new function node with the same head. Returns `NULL`
(unevaluated) for atoms or empty expressions.
