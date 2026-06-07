---
source: src/core.c
---
`builtin_compose_list` (`src/core.c`) takes `{f1,...,fn}` and `x` and builds the length-`n+1` list `{x, f1[x], f2[f1[x]], ...}` by constructing each symbolic application `fi[prev]`; the outer evaluator then reduces those applications to fixed point. Returns `NULL` if the first argument is not a `List`.
