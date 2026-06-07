# SubresultantPolynomials

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SubresultantPolynomials[poly1, poly2, var]
    gives the list of subresultant polynomials {S_0, ..., S_m} of
    poly1 and poly2 with respect to var, where m = Exponent[poly2,
    var].  The list has length m + 1, its first element is
    Resultant[poly1, poly2, var], and the coefficient of var^j in S_j
    is the j-th principal subresultant coefficient.  Requires
    Exponent[poly1, var] >= Exponent[poly2, var] and exact
    coefficients.  Computed by a subresultant polynomial-remainder
    sequence, with a determinant-polynomial fallback for algebraic
    coefficients.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SubresultantPolynomials[(x - 1)^2 (x - 2) (x - 3), (x - 1) (x - 4)^2, x]
Out[1]= {0, -36 + 36 x, 38 - 49 x + 11 x^2, -16 + 24 x - 9 x^2 + x^3}

In[2]:= SubresultantPolynomials[a x^3 + b x^2 + c x + d, 3 a x^2 + b x + c, x]
Out[2]= {4 a^2 c^3 + 2 a b^3 d - 18 a^2 b c d + 27 a^3 d^2, -2 a b c + 9 a^2 d - 2 a b^2 x + 6 a^2 c x, c + b x + 3 a x^2}

In[3]:= SubresultantPolynomials[2 x^7 + 3 x^3 + 5 x - 1, 7 x^6 + 8 x - 9, x]
Out[3]= {-183782157189, -761749829 + 3208696817 x, -3143546 + 11222638 x + 3838135 x^2, -21609 + 163611 x - 49392 x^2 + 64827 x^3, 0, -49 + 371 x - 112 x^2 + 147 x^3, -9 + 8 x + 7 x^6}

In[4]:= First[%] - Resultant[2 x^7 + 3 x^3 + 5 x - 1, 7 x^6 + 8 x - 9, x]
Out[4]= 0
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
