---
source: src/part.c
---
**Algorithm.** `builtin_part` extracts elements by index path, delegating to the recursive
`expr_part(expr, indices, nindices)`. Each index level may be: a positive or negative integer
(`-k` resolves to `len + k + 1`); `0`, which extracts the head and is allowed even on atoms; a
`List` of indices (extract several, returning a list); `All`; or a `Span` (`i;;j;;k`) built by
the parser from `;;` syntax, which `expr_part` resolves into an explicit element range with the
given start/end/step (negative endpoints wrap, `UpTo`/`All` endpoints clamp). Index paths apply
left to right, descending one structural level per index. Out-of-range or non-integer indices on
atoms yield `NULL` (unevaluated).
