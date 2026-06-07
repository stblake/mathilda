---
source: src/poly/poly.c
---
**Algorithm.** `builtin_variables` walks the expression with `collect_variables`, gathering the
distinct symbols that occur as polynomial generators (bare symbols and non-numeric bases,
excluding numeric atoms and known constants), then sorts the collected `Expr*` array with
`qsort` under `compare_expr_ptrs` (the canonical `expr_compare` order) and wraps the result in a
`List`. The output is the deduplicated, canonically-ordered list of variables on which the
input is treated as a polynomial/rational expression.
