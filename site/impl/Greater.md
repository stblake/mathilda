---
source: src/comparisons.c
---
`builtin_greater` delegates to the shared `evaluate_inequality(res, 1, 1)`, which walks adjacent argument pairs and calls `compare_numeric` (exact GMP for integer-like, exact cross-multiplied `long double` for rationals, `2^-46` relative tolerance for inexact reals; returns -1/0/+1). A pair is accepted only when the sign is `+1` (strictly greater). All pairs strictly decreasing → `True`; any `0`/`-1` pair → `False`; any non-comparable (symbolic) pair → NULL (chain left unevaluated). This is the chained `a > b > c` semantics. Shares its implementation with `Less`/`LessEqual`/`GreaterEqual` via different accepted-sign arguments.
