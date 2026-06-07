---
source: src/boolean.c
---
`builtin_boole` is a one-argument map from boolean symbols to integers: it returns `1` for the interned `True`, `0` for `False`, and `NULL` (stays unevaluated) for anything else. It is registered `ATTR_LISTABLE | ATTR_PROTECTED`, so the evaluator threads it over `List` arguments automatically and the handler only needs the scalar case.
