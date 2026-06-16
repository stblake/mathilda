### Worked examples

```mathematica
In[1]:= PolynomialQuotientRemainder[x^2 - 1, x - 1, x]
Out[1]= {1 + x, 0}
```

```mathematica
In[1]:= PolynomialQuotientRemainder[x^5 + x + 1, x^2 + 1, x]
Out[1]= {-x + x^3, 1 + 2 x}
```

```mathematica
In[1]:= {q, r} = PolynomialQuotientRemainder[x^5 + x + 1, x^2 + 1, x];
In[2]:= Expand[q (x^2 + 1) + r]
Out[2]= 1 + x + x^5
```

```mathematica
In[1]:= PolynomialQuotientRemainder[x^4 - 2, x^2 - Sqrt[2], x, Extension -> Sqrt[2]]
Out[1]= {Sqrt[2] + x^2, 0}
```

### Notes

`PolynomialQuotientRemainder[p, q, x]` performs a single long division and
returns `{quotient, remainder}` together, satisfying `p == quotient*q +
remainder` with `deg(remainder) < deg(q)`. The third example reconstructs the
dividend `x^5 + x + 1` from the returned pair, verifying the division identity.
The `Extension -> alpha` option carries out the division over `Q(alpha)[x]`;
over `Q(Sqrt[2])` the polynomial `x^4 - 2 = (x^2 - Sqrt[2])(x^2 + Sqrt[2])`
divides exactly, giving a zero remainder. This is the combined form of
`PolynomialQuotient` and `PolynomialRemainder`, computed in one pass.
