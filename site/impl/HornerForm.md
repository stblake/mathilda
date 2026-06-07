---
source: src/poly/poly.c
references:
  - "W. G. Horner, \"A new method of solving numerical equations of all orders, by continuous approximation\", Phil. Trans. R. Soc. 1819."
---
**Algorithm.** `builtin_hornerform` (in `src/poly/poly.c`) rewrites a polynomial in nested multiply-add form. It first splits a rational input into numerator/denominator (scanning a top-level `Times`/`Power` for negative-exponent factors), then applies the recursive worker `horner_form_rec` to the numerator. That worker `expr_expand`s, verifies the expression is a polynomial in the lead variable via `PolynomialQ`, obtains its dense coefficients with `CoefficientList`, and folds from the highest-degree coefficient downward using the Horner recurrence `H ← c_i + v·H` (with zero-coefficient short-circuits). For multiple variables it recurses on each coefficient with the remaining variable list, producing a nested Horner form.

**Data structures.** Reuses the `CoefficientList` array as the coefficient sequence; nested `Times`/`Plus` nodes (built via `internal_times`/`internal_plus`) form the output.
