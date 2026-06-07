---
references:
  - "G. E. Collins, \"Subresultants and Reduced Polynomial Remainder Sequences\", JACM 14(1), 1967."
source: src/poly/poly.c
---
**Algorithm.** `builtin_polynomiallcm` mirrors `PolynomialGCD`: it strips an optional
`Extension -> α`/`Automatic` (routing to `polynomiallcm_with_extension` /
`qa_polynomiallcm_with_tower*` over `Q(α)`), force-rationalises and re-numericalises inexact
coefficients, and otherwise pre-decomposes each input with `decompose_to_bp` to handle numeric
content and shared non-numeric factors. The polynomial parts are combined pairwise via the
identity `lcm(a, b) = a·b / gcd(a, b)`, where the GCD comes from `poly_gcd_internal` (the
recursive multivariate subresultant PRS — see `PolynomialGCD`) and the exact division is done
by `exact_poly_div`/`Cancel`. Multiple arguments fold left-to-right, accumulating the running
LCM. When the exact quotient by the GCD can't be certified the code falls back to the plain
product `a·b`.

**Data structures.** `Expr` trees throughout; `BPList` for the base/power decomposition;
`QAUPoly`/`QATower` on the algebraic-extension path.
