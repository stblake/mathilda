---
source: src/comparisons.c
---
`builtin_sameq` tests purely structural identity. With fewer than two arguments it returns `True` by convention; otherwise it compares every argument against the first with `expr_eq` and returns `True` only if all are structurally equal, `False` otherwise. Unlike `Equal`, there is no numeric coercion — `1 === 1.` is `False` — and it never returns `NULL` (the result is always a definite boolean).
