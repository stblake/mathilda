---
source: src/poly/poly.c
---
**Algorithm.** `builtin_polynomialremainder` is the companion to `PolynomialQuotient`: after
stripping an optional `Extension -> α`/`Automatic` option (routing to
`polynomialdivrem_with_extension` over `Q(α)`), it calls the shared
`poly_div_rem(p, q, x, &Q, &R)` and returns the remainder `R`, discarding the quotient. The
division is field long division in `x` (see `PolynomialQuotient`): subtract
`(lc(R)/lc(q))·x^(deg R − deg q)·q` from the running remainder until `deg R < deg q`, with a
fast exact-integer-division path for pure integer/bigint leading coefficients. The
`PolynomialQuotientRemainder` builtin exposes both outputs `{Q, R}` from a single call.

**Data structures.** Expanded `Expr` polynomial trees; coefficients extracted via
`get_coeff_expanded`.
