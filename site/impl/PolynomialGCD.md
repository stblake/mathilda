---
references:
  - "W. S. Brown, \"On Euclid's Algorithm and the Computation of Polynomial Greatest Common Divisors\", JACM 18(4), 1971."
  - "G. E. Collins, \"Subresultants and Reduced Polynomial Remainder Sequences\", JACM 14(1), 1967."
source: src/poly/poly.c
---
**Algorithm.** `builtin_polynomialgcd` first strips an optional `Extension -> α` (or
`Extension -> Automatic`) and, when present, routes through the algebraic-number machinery
(`polynomialgcd_with_extension`, which lifts each input to a `QAUPoly` over `Q(α)` and folds
them with `qaupoly_gcd`; `qa_polynomialgcd_with_tower*` for multi-generator towers). Inexact
(floating) coefficients are force-rationalised, run through the exact algorithm, then
numericalised (`internal_rationalize_then_numericalize`).

The core path pre-processes each input with `decompose_to_bp` into a base/power list, peeling
off the integer content (numeric GCD of literal coefficients, including the integer content of
`Plus`-headed factors so it isn't double-counted) and any non-numeric factors common to every
argument. The remaining symbolic GCD is computed by `poly_gcd_internal`, which implements the
**recursive multivariate subresultant PRS**: it treats the last variable as main, splits each
operand into content (GCD of its coefficients, computed recursively in one fewer variable) and
primitive part, then reduces the primitive parts with `pseudo_rem` (a pseudo-remainder that
stays inside the coefficient ring, avoiding rationals) until the remainder is zero. The base
case (zero variables) is integer GCD via `my_number_gcd`. The result is `content_GCD ×
primitive_GCD`, normalised to a positive leading coefficient and expanded. Multi-argument GCD
folds left-to-right. A size budget (`max(input_size, 2000)` leaves) and a 50-iteration cap
guard against coefficient explosion over multi-radical rings; on overflow it conservatively
returns just the content GCD (always a valid divisor).

**Data structures.** Inputs are `Expr` trees; `BPList` holds the base/power decomposition;
`QAUPoly`/`QAExt`/`QATower` carry the algebraic-extension representation. Coefficients are
ordinary `Expr` subtrees, so coefficient GCDs recurse through the same machinery.
