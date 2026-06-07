---
source: src/part.c
---
`builtin_first` (in `src/part.c`) takes a single argument and returns a deep copy of its first element (`args[0]`). It returns `NULL` (unevaluated) when the argument is atomic or has no elements.
