---
source: src/part.c
---
`builtin_insert` (in `src/part.c`) handles the 3-arg form `Insert[expr, elem, pos]` by delegating to the helper `expr_insert`, which inserts `elem` before position `pos` (negative indices count from the end, and a position may be a list path for nested insertion), deep-copying the surrounding structure and preserving the original head.
