---
source: src/list.c
---
`builtin_square_matrix_q` is a pure shape test: it returns `True` iff the argument is a non-empty `List` of equal-length `List`s with row count equal to column count and no entry that is itself a `List` (rejecting ragged, rectangular, and higher-rank tensors). No element predicate is consulted, so `{{x}}` is square for any `x`. Exactly one argument is accepted; any other count emits `SquareMatrixQ::argx` and leaves the call unevaluated.
