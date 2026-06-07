---
source: src/linalg/det.c
---
**Algorithm.** `builtin_det` validates that the argument is a non-empty square rank-2 tensor (via `get_tensor_dims`; otherwise it emits `Det::matsq` and returns `NULL`), flattens it row-major into an `Expr**`, and computes the determinant by full **Laplace cofactor expansion** along the first row, recursively (`laplace_det`). Each cofactor term is built as `Times[±1, element, minor]` and accumulated with `Plus`, every product/sum being reduced through `eval_and_free`, so cancellation and symbolic simplification happen as the expansion unwinds. This keeps results exact and symbolic for integer, rational, and symbolic matrices.

**Data structures.** The matrix is a flat `Expr**` of `n*n` element pointers; recursion carries an explicit `int* cols` index set and a fixed `row` cursor, deleting one column per level rather than copying submatrices. `laplace_det` is exported via `linalg.h` and reused by `Cross`.

**Complexity / limits.** Cofactor expansion is `O(n!)`, so it is only practical for small `n`. There is no fraction-free Bareiss or LU fast path in this handler — the larger fraction-free Gauss-Jordan machinery lives in `inv.c`/`linsolve.c` for inversion and solving, not in `Det`.
