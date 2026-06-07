---
source: src/comparisons.c
---
`builtin_unsameq` is the negation of structural identity over all argument pairs. With fewer than two arguments it returns `True`; otherwise it checks every pair `(i,j)` with `expr_eq` and returns `False` as soon as any two are structurally equal, else `True`. Like `SameQ` it uses no numeric coercion and always yields a definite boolean.
