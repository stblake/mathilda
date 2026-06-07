---
source: src/linalg/nullspace.c
---
**Algorithm.** `builtin_nullspace` computes a basis for `{v : m.v == 0}` by reduction to RREF rather than by an orthogonal-completion (QR/SVD) route. `nullspace_core` calls `RowReduce[m, Method -> ...]` through the evaluator (`call_rowreduce`, forwarding the optional `Method` option), flattens the RREF, and locates each row's pivot column as its leftmost structurally non-zero entry (`is_zero_poly`). For every *free* (non-pivot) column `f`, iterated rightmost-to-leftmost to match Mathematica's ordering, it builds a length-`cols` basis vector with `v[f] = 1`, `v[p] = -RREF[row_of_p, f]` for each pivot column `p`, and `0` elsewhere.

**Data structures.** A flat `Expr**` of the RREF plus an `int* pivot_for_col` map; basis vectors are accumulated into a growable `Expr**`. For exact rational input each vector is scaled by the LCM of its entries' integer denominators (`clear_int_denominators`, via GMP `mpz_lcm` and the `Denominator` builtin) so an integer-valued nullspace comes out integer-valued; symbolic/inexact entries are left in natural form. Full-column-rank input returns `List[]`.

**Limits.** Inherits whatever the `RowReduce` dispatcher (default `DivisionFreeRowReduction`, Bareiss-like fraction-free Gauss-Jordan) can reduce; correctness for symbolic matrices depends on `is_zero_poly`/`Together` detecting cancellations. Bad `Method` values emit `NullSpace::method`; non-rectangular input emits `NullSpace::matrix`.
