### Worked examples

The subresultant polynomials `{S_0, ..., S_m}` of two quadratics; `S_0` is the
resultant and the top entry is the leading polynomial:

```mathematica
In[1]:= SubresultantPolynomials[x^2 - 1, x^2 - 4, x]
Out[1]= {9, -3, 1}
```

A cubic/quadratic pair, where `S_2` is itself a genuine polynomial in `x`:

```mathematica
In[1]:= SubresultantPolynomials[x^3 + x + 1, x^2 + 1, x]
Out[1]= {1, 1, 1 + x^2}
```

`S_0` agrees with `Resultant` as expected:

```mathematica
In[1]:= Resultant[x^2 - 1, x^2 - 4, x]
Out[1]= 9
```

### Notes

`SubresultantPolynomials[p1, p2, x]` returns the full chain of subresultant
*polynomials* `{S_0, ..., S_m}` with `m = Exponent[p2, x]`, so the list has
length `m + 1`. Its first element is `Resultant[p1, p2, x]`, and the coefficient
of `x^j` in `S_j` is the `j`-th principal subresultant coefficient (the value
returned by `Subresultants`). It requires `Exponent[p1, x] >= Exponent[p2, x]`
and exact coefficients, and is computed by a subresultant polynomial-remainder
sequence with a determinant-polynomial fallback for algebraic coefficients.
