---
source: src/linalg/hankelmat.c
---
**Algorithm.** `builtin_hankelmatrix` builds a matrix constant along antidiagonals, where entry `(i,j)` depends only on `s = i+j-1`. Three forms are handled: `HankelMatrix[n]` (square, antidiagonal index `s` for `s ≤ n` else `0`, the integer form); `HankelMatrix[{c1,…,cm}]` (square, first column `c`, zeros below the antidiagonal); and `HankelMatrix[{c…}, {r…}]` (`m×n`, first column `c` and last row `r`, with `(i,j) = c_s` when `s ≤ m` and `r_{s-m}` otherwise). The shared corner `c_m` must equal `r_1`; if not, it warns via `HankelMatrix::crs` and uses the column element. Zero arguments emit `HankelMatrix::argb`; any other shape (non-list, empty list) is left unevaluated.

**Data structures.** A `List` of `List`s built by `hk_build`; source entries are deep-copied (`expr_copy`) so symbolic/complex/exact/inexact entries pass through verbatim — arbitrary precision comes from the entries themselves. Complexity `O(mn)`.
