---
source: src/poly/poly.c
---
**Algorithm.** `builtin_coefficientlist` (in `src/poly/poly.c`) returns the dense coefficient array of a polynomial. It `expr_expand`s the input, computes each variable's degree with `get_degree_poly`, then the recursive worker `coeff_list_rec` builds a nested `List` whose shape mirrors the variable order. At each level it pulls all coefficients `c_0..c_d` of the current variable — preferring the bulk extractor `get_all_coeffs_expanded` (single pass over the expanded form) and falling back to `get_coeff_expanded` per degree — and recurses on each coefficient for the next variable.

**Data structures.** A `int* max_degrees` array sizes each axis; coefficients are produced as `Expr*` and assembled into nested `List` nodes. The bulk path avoids the naive (degree+1)-passes-per-level cost.
