---
references:
  - "K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992) — Ch. 2, polynomial division."
source: src/poly/poly.c
---
**Algorithm.** `PolynomialQuotientRemainder[p, q, x]` returns the pair
`{quotient, remainder}` from Euclidean division of `p` by `q` in `x`, computing
both halves in one shot. `builtin_polynomialquotientremainder`
(src/poly/poly.c) calls the shared `poly_div_rem` helper, which expands both
operands, reads the leading coefficient of `q`, and runs the classical
long-division loop: at each step it forms the quotient term
`lc(R)/lc(q) · x^(deg R − deg q)`, subtracts `term·q` from the running remainder,
and repeats until `deg R < deg q`. Constant divisors are short-cut to
`{p/q, 0}`. The quotient is `Expand`-ed before being returned in a `List`.

An optional `Extension -> α` (or `Extension -> Automatic`, which autodetects the
algebraic generators of `p` and `q` via `extension_autodetect_args`) re-runs the
division inside an algebraic number field tower (`QATower`,
`polynomialdivrem_with_extension`), falling back to the plain path on lift
failure.

**Data structures.** `Expr*` polynomial trees; coefficients are exact
(`EXPR_INTEGER`/`EXPR_BIGINT`/`Rational`), with an integer/bigint fast path
(`mpz_tdiv_qr`) that avoids the `Together`/`Cancel` denominator unification when
coefficients stay integral. Extension arithmetic uses the `QATower` algebraic
number representation.

**Complexity / limits.** `O(deg p · deg q)` coefficient operations for the
machine path. The divisor must be nonzero (returns NULL otherwise) and the
variable a symbol. Multivariate inputs are handled coefficient-wise in `x`.
