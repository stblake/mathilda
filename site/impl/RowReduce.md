---
references:
  - "E. H. Bareiss, \"Sylvester's Identity and Multistep Integer-Preserving Gaussian Elimination\", Math. Comp. 22 (1968)."
source: src/linalg/linsolve.c
---
**Algorithm.** `builtin_rowreduce` reduces `m` to reduced row-echelon form, dispatching on an optional `Method` option (`matsol_parse_method_option`):

- **`DivisionFreeRowReduction`** (the default / `Automatic`): `rowreduce_divfree`, a Bareiss-like fraction-free Gauss-Jordan. Per pivot column it picks the first structurally non-zero pivot (`is_zero_poly`), eliminates every other row with the Bareiss update `M[i,j] ← (pivot·M[i,j] − M[i,c]·M[r,j]) / P` where `P` is the previous pivot and the division is exact (`exact_div_wrapper`), keeping all intermediates polynomial/integer with no GCD blow-up. A final pass scales each pivot row to a leading 1, reducing each entry by its `PolynomialGCD` with the leading coefficient and `expand`-ing.
- **`OneStepRowReduction`**: `rowreduce_onestep`, textbook Gauss-Jordan — normalise the pivot row by dividing by the pivot (`matsol_div_entry`, which prefers exact division and canonicalises with `Together`), then subtract multiples from every other row. Result is already RREF.
- **`CofactorExpansion`**: `rowreduce_cofactor` returns identity-if-invertible for a non-singular square matrix, otherwise warns (`RowReduce::cofnsq`) and falls back to the fraction-free path.

**Data structures.** A flat row-major `Expr**` of size `m·n`; all arithmetic (`Times`, `Plus`, `Power`, `PolynomialGCD`, `Together`) is driven through the evaluator via `eval_and_free`, so integer, rational, complex, and free-symbolic matrices share one code path. Pivot detection is purely structural (`is_zero_poly`), so a `Real` cancellation of magnitude ~1e-18 is not treated as zero. Bad `Method` values emit `RowReduce::method` and leave the call unevaluated.
