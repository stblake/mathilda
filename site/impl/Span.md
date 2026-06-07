---
source: src/part.c
---
**Algorithm.** `Span` is not a standalone builtin — it is the head the parser produces for the
`i;;j` and `i;;j;;k` syntax (`OP_SPAN`, precedence 290, in `parse.c`). A `Span[start, end,
step]` expression is interpreted only inside `Part`/`Extract`: `expr_part` (in `part.c`)
resolves it against a concrete list, mapping negative endpoints to `len + k + 1`, clamping
`UpTo`/`All` endpoints, computing the element count from the step, and emitting the selected
elements. On its own a bare `Span` stays unevaluated.
