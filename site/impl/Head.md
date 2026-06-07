---
source: src/part.c
---
`builtin_head` (in `src/part.c`) returns the head of its argument via `expr_head` — `f` for `f[...]`, and the type symbol (`Integer`, `Symbol`, `List`, etc.) for atoms. The 2-arg form `Head[expr, h]` wraps the result as `h[Head[expr]]`, leaving the outer application for the evaluator to reduce.
