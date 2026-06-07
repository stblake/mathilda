---
source: src/part.c
---
`builtin_last` (in `src/part.c`) takes a single argument and returns a deep copy of its final element (`args[arg_count - 1]`). It returns `NULL` (unevaluated) when the argument is atomic or empty.
