---
source: src/core.c
---
`builtin_append` (in `src/core.c`) requires a 2-arg call `Append[expr, elem]` whose first argument is an `EXPR_FUNCTION`. It allocates a fresh argument array one slot longer than `expr`, deep-copies every existing argument plus `elem` into it, and rebuilds a new `EXPR_FUNCTION` with the same head. Works on any head, not just `List`. Returns `NULL` (unevaluated) when the first argument is atomic.
