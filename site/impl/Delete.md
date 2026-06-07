---
source: src/part.c
---
`builtin_delete` (in `src/part.c`) drives the recursive helper `delete_path`, which walks an integer position (or position path) into the expression tree and removes the targeted element by rebuilding the enclosing function with all arguments except that index. Negative indices count from the end, position `0` targets the head (replaced by a `Sequence[...]` of the remaining parts), and out-of-range indices leave the structure unchanged.
