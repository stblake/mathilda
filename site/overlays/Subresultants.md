### Worked examples

The principal subresultant coefficients of two coprime quadratics; the leading
entry is the resultant:

```mathematica
In[1]:= Subresultants[x^2 - 1, x^2 - 4, x]
Out[1]= {9, 0, 1}
```

That first entry matches `Resultant` exactly:

```mathematica
In[1]:= Resultant[x^2 - 1, x^2 - 4, x]
Out[1]= 9
```

A coprime cubic/quadratic pair:

```mathematica
In[1]:= Subresultants[x^3 + x + 1, x^2 + 1, x]
Out[1]= {1, 0, 1}
```

When the polynomials share a root the resultant vanishes — the leading entries
are exactly zero, signalling the common factor `(x - 1)`:

```mathematica
In[1]:= Subresultants[x^2 - 1, (x - 1)^2, x]
Out[1]= {0, -2, 1}
```

### Notes

`Subresultants[p1, p2, x]` returns the list of principal subresultant
coefficients with respect to `x`. Its length is
`Min[Exponent[p1, x], Exponent[p2, x]] + 1`, its first element equals
`Resultant[p1, p2, x]`, and its first `k` entries vanish exactly when the
polynomials share `k` roots (counted with multiplicity) — making it a robust
GCD-degree detector that needs no factorization. It is computed from a
subresultant polynomial-remainder sequence (with a Sylvester-minor determinant
fallback for algebraic coefficients).
