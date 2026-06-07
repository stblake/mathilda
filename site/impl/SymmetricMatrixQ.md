---
source: src/list.c
---
`builtin_symmetric_matrix_q` first applies the same square-matrix shape gate as `SquareMatrixQ`, then walks the strict upper triangle checking `m[i,j] == m[j,i]`. The comparison defaults to structural `expr_eq`, but a `SameTest -> f` option uses `symmetric_pair_sametest` and a `Tolerance -> t` option uses `symmetric_pair_tolerance`. Returns `False` on any shape rejection or mismatch; unrecognised options leave the call unevaluated.
