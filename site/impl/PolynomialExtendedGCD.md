---
references:
  - "Extended Euclidean algorithm over a polynomial ring; Bézout's identity."
source: src/poly/poly.c
---
**Algorithm.** `builtin_polynomialextendedgcd` runs the standard **extended Euclidean
iteration** on the two univariate inputs `A`, `B` in `x`, maintaining the cofactor triples and
updating each pass by `r_{i+1} = r_{i-1} − q_i r_i`, `s_{i+1} = s_{i-1} − q_i s_i`, `t_{i+1} =
t_{i-1} − q_i t_i`, where the quotient `q_i` comes from polynomial division. On exit `r_0` is
the GCD and `(s_0, t_0)` are the Bézout cofactors satisfying `s·A + t·B = gcd`. The GCD is
finally normalised to be monic in `x` (dividing the triple through by the leading coefficient),
matching Mathematica's `{g, {s, t}}` result shape. An optional fourth argument `Modulus -> p`
switches the coefficient arithmetic to `Z/pZ`, using `mod_inverse_int_poly` (itself the
extended Euclidean algorithm on integers) to invert leading coefficients.

**Data structures.** Operands and cofactors are ordinary `Expr` polynomial trees in `x`;
division/remainder reuse the field-based `poly_div_rem` long-division routine.
