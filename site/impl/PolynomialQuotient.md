---
source: src/poly/poly.c
---
**Algorithm.** `builtin_polynomialquotient` strips an optional `Extension -> α`/`Automatic`
(routing to `polynomialdivrem_with_extension` over `Q(α)` when applicable), then calls the
shared `poly_div_rem(p, q, x, &Q, &R)` and returns the expanded quotient `Q`, discarding `R`.
`poly_div_rem` is textbook **field long division** in `x`: it expands both operands, repeatedly
forms the next quotient term `(lc(R)/lc(q)) · x^(deg R − deg q)`, subtracts `term·q` from the
running remainder `R`, and stops when `deg R < deg q`. A fast path detects exact integer/bigint
leading-coefficient divisions (via `mpz_tdiv_qr`) so the subtraction step never needlessly
introduces rationals; otherwise quotient coefficients are formed symbolically. Coefficients are
extracted with `get_coeff_expanded` against the already-expanded divisor.

**Data structures.** `Expr` polynomial trees in expanded form; quotient and remainder are
returned through `out_Q`/`out_R` pointers.
