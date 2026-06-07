---
source: src/linalg/matrank.c
---
**Algorithm.** `builtin_matrixrank` returns the rank as the number of pivots, choosing between two paths. The **numerical path** is taken when the user supplies a finite `Tolerance` or the matrix has any inexact (`Real`/`MPFR`) leaf: every entry is coerced to a `cplx_t` `{re, im}` struct via `entry_to_cplx` (handling Integer, BigInt, Real, MPFR, `Rational`, `Complex`, `I`, and a `N[expr]` fallback for `Pi`, `Sqrt[2]`, etc.), then `gauss_rank_cplx` runs partial-pivot Gaussian forward elimination over `double`-complex and counts accepted pivots — a column is skipped when its largest sub-pivot magnitude is `<= tol`. The default tolerance for inexact input is `max(rows,cols) · DBL_EPSILON · max(|entries|)` (the standard rank-by-SVD surrogate, with `max(|entries|)` substituted for the top singular value); for exact input the default is `0`.

The **exact path** (every leaf exact, no `Tolerance`, or the numeric coercion failed because of symbolic entries) calls `RowReduce[m, Method -> ...]` through the evaluator (`call_rowreduce`) and counts RREF rows that are not all structurally zero (`count_nonzero_rows`, using `is_zero_poly`). An optional `Method` option (`DivisionFreeRowReduction`, `OneStepRowReduction`, `CofactorExpansion`, `Automatic`) is forwarded to `RowReduce`.

**Data structures / limits.** Numerical path: a flat `cplx_t` array of size `rows*cols` with hand-rolled complex add/sub/mul/div. Exact path: defers entirely to the `RowReduce` dispatcher. Bad options emit `MatrixRank::opt`; non-rectangular input emits `MatrixRank::matrix`.
