---
source: src/list.c
---
`builtin_diagonal_matrix_q` tests whether a matrix has nonzero entries only on the `k`-th diagonal (default `k = 0`). It accepts an optional integer `k` at position 2 and a `Tolerance` option; an empty or malformed argument list yields a `DiagonalMatrixQ::argt`/`::nonopt` diagnostic, and missing args return `False`. After validating that the matrix is a rectangular `List` of `List`s, it returns `True`/`False` according to whether every off-`k`-diagonal entry is structurally (or within `Tolerance`) zero.
