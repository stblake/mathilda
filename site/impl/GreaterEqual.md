---
source: src/comparisons.c
---
`builtin_greaterequal` delegates to `evaluate_inequality(res, 1, 0)`: it walks adjacent argument pairs through `compare_numeric` (exact GMP for integer-like operands, exact cross-multiplied `long double` for rationals, `2^-46` relative tolerance for inexact reals) and accepts a pair when the sign is in `{1, 0}` (greater-than-or-equal). All pairs satisfying ≥ → `True`; any pair comparing `-1` → `False`; any non-comparable symbolic pair → NULL (left unevaluated). Implements chained `a >= b >= c`. Same implementation as `Less`/`Greater`/`LessEqual` with a different accepted-sign set.
