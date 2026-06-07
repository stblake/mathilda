---
source: src/cond.c
---
`builtin_trueq` is a one-argument predicate (no Hold attributes, so its argument is already evaluated). It returns the symbol `True` only when the argument is exactly the interned symbol `True` (pointer equality on `SYM_True`), and `False` for everything else — so unlike a bare condition it never stays symbolic. A non-unary call returns `NULL`.
