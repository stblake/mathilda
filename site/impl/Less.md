---
source: src/comparisons.c
---
`builtin_less` delegates to the shared `evaluate_inequality(res, -1, -1)`, which walks adjacent argument pairs and calls `compare_numeric` on each. `compare_numeric` returns -1/0/+1 using exact GMP comparison for integer-like operands, exact cross-multiplied `long double` comparison for rationals, and a `2^-46` relative-tolerance comparison for inexact reals. `Less` accepts a pair only when the sign is in `{-1, -1}` (strictly less). If every pair is strictly increasing → `True`; any pair that compares `0` or `+1` → `False`; any pair that is not numerically comparable → NULL (the whole chain stays unevaluated). This implements the chained semantics of `a < b < c`. `Greater`/`LessEqual`/`GreaterEqual` are the same function with different accepted sign sets.
