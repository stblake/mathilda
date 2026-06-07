---
source: src/linalg/toeplitzmat.c
---
**Algorithm.** `builtin_toeplitzmatrix` constructs a matrix constant along diagonals. Three forms dispatch on argument shape: `ToeplitzMatrix[n]` builds the symmetric `n × n` integer matrix with entry `(i,j) = |i−j| + 1`; `ToeplitzMatrix[{c}]` builds the symmetric matrix whose first column and row are the list `c`; `ToeplitzMatrix[{c}, {r}]` builds the `m × n` matrix with `c` down the first column and `r` across the first row. The builder `tz_build` sets entry `(i,j)` to `cvals[i−j]` when `i >= j` and `rvals[j−i]` otherwise, deep-copying source entries verbatim so symbolic/complex/exact/inexact entries flow through unchanged.

**Limits.** The shared corner reads `c_1` on the diagonal; if `c_1 != r_1` it warns `ToeplitzMatrix::crs` and uses the column element. Zero arguments emit `ToeplitzMatrix::argb`; empty lists or any other shape leave the call unevaluated.
