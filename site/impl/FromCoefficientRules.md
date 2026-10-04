---
source: src/poly/monomials.c
---
**Algorithm.** `builtin_fromcoefficientrules` reconstructs a polynomial from the
sparse form `CoefficientRules` produces — it is the inverse of that head.
`FromCoefficientRules[{expvec -> coeff, ...}, {x1, ..., xk}]` reads the variable
list (a `List` or a single bare variable, giving `k`), then walks each rule: the
left side must be a length-`k` `List` of `Integer` exponents and the right side
is the coefficient. For each rule it assembles the monomial
`coeff * x1^e1 * ... * xk^ek` with `internal_power`/`internal_times`, skipping a
variable whose exponent is `0` and dropping the `Power` wrapper when the exponent
is `1`. The monomials are summed with `internal_plus`; an empty rule list
reconstructs the `Integer` `0`.

**Data structures.** No intermediate monomial table — each rule is turned
straight into a `Times` node in a growable `Expr**` `terms` buffer, then folded
into a single `Plus`. A malformed rule (wrong head, a non-`List` or wrong-length
exponent vector, or a non-`Integer` exponent) frees the partial terms and leaves
the call unevaluated (`NULL`).

**Complexity / limits.** Linear in the number of rules times `k`; the real cost
is the evaluator's canonicalisation of the assembled `Plus` (collecting like
terms). A symbolic structural head — exact `Expr` arithmetic, so no
packed/NDArray or `Compile[]` path. `Protected`.
