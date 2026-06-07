---
source: src/list.c
---
`builtin_upper_triangular_matrix_q` validates that the argument is a non-empty rectangular `List` of equal-length `List`s (no deeper nesting), then returns `True` iff every entry strictly below the `k`-th diagonal (column−row `< k`, default `k = 0`) is zero. The optional second Integer argument selects the diagonal `k`; a `Tolerance -> t` option relaxes the zero test (otherwise a structural zero check is used). Bad arguments/options emit `UpperTriangularMatrixQ::nonopt` / `::argt`; shape rejections return `False`.
