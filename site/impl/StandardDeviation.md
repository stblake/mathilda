---
source: src/stats.c
---
`builtin_standard_deviation` is essentially `Sqrt[Variance[data]]`. It reduces matrices column-wise via `apply_columnwise`. For an all-real numeric vector (`n > 1`) it evaluates `Variance[data]` and, if that returns an `EXPR_REAL`, returns `expr_new_real(sqrt(...))` directly. Otherwise it evaluates `Variance[data]` and raises it to the `1/2` power via a `Power[var, Rational[1,2]]` node, letting the evaluator produce an exact or symbolic radical. `ATTR_PROTECTED`. Inherits Variance's sample (`n-1`) convention.
