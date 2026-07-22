# PolynomialExtendedGCD

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PolynomialExtendedGCD[poly1, poly2, x] gives the extended GCD of poly1 and poly2 treated as univariate polynomials in x.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PolynomialExtendedGCD[2x^5-2x, (x^2-1)^2, x]
Out[1]= {-1 + x^2, {1/4 x, 1/2 (-2 - x^2)}}

In[2]:= PolynomialExtendedGCD[a (x+b)^2, (x+a)(x+b), x]
Out[2]= {b + x, {-1/(a^2 - a b), 1/(a - b)}}
```

## Implementation notes

**Algorithm.** `builtin_polynomialextendedgcd` runs the standard **extended Euclidean
iteration** on the two univariate inputs `A`, `B` in `x`, maintaining the cofactor triples and
updating each pass by `r_{i+1} = r_{i-1} − q_i r_i`, `s_{i+1} = s_{i-1} − q_i s_i`, `t_{i+1} =
t_{i-1} − q_i t_i`, where the quotient `q_i` comes from polynomial division. On exit `r_0` is
the GCD and `(s_0, t_0)` are the Bézout cofactors satisfying `s·A + t·B = gcd`. The GCD is
finally normalised to be monic in `x` (dividing the triple through by the leading coefficient),
matching the `{g, {s, t}}` result shape. An optional fourth argument `Modulus -> p`
switches the coefficient arithmetic to `Z/pZ`, using `mod_inverse_int_poly` (itself the
extended Euclidean algorithm on integers) to invert leading coefficients.

**Data structures.** Operands and cofactors are ordinary `Expr` polynomial trees in `x`;
division/remainder reuse the field-based `poly_div_rem` long-division routine.

- `Protected`.
- Returns `{d, {a, b}}` such that $a \cdot poly1 + b \cdot poly2 = d$.
- $d$ is the GCD, normalized to be monic.
- Efficiently handles termination when a constant remainder is reached.
- Optimized for cases where the divisor is a constant.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Extended Euclidean algorithm over a polynomial ring; Bézout's identity.
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= PolynomialExtendedGCD[x^2 - 1, x^3 - 1, x]
Out[1]= {-1 + x, {-x, 1}}
```

```mathematica
In[1]:= PolynomialExtendedGCD[x^4 + x^3 + x^2 + x + 1, x^2 + 1, x]
Out[1]= {1, {1, -x - x^2}}
```

```mathematica
In[1]:= PolynomialExtendedGCD[x^7 - 1, x^5 - 1, x]
Out[1]= {-1 + x, {-x - x^3, 1 + x^3 + x^5}}
```

### Notes

`PolynomialExtendedGCD[a, b, x]` returns `{g, {s, t}}` where `g` is the
polynomial GCD and `s`, `t` are the Bezout cofactors satisfying
`s a + t b == g`. The first example certifies `(-x)(x^2-1) + 1*(x^3-1) = x - 1`,
the GCD. When the inputs are coprime the GCD is `1` and the cofactors give a
constructive proof of coprimality — exactly the data needed to invert one
polynomial modulo another (e.g. building inverses in `Q[x]/(b)`). For the two
cyclotomic-flavoured inputs `x^7-1` and `x^5-1`, whose only common root is the
trivial seventh-and-fifth root `x = 1`, the GCD is the linear factor `x - 1`.
