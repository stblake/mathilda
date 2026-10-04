### Worked examples

```mathematica
In[1]:= PolynomialReduce[x^2 + y^2, {x - y}, {x, y}]  (* {{quotient}, remainder} *)
```

```mathematica
In[2]:= PolynomialReduce[x^2 y + x y^2 + y^2, {x y - 1, y^2 - 1}, {x, y}]  (* several divisors *)
```

```mathematica
In[3]:= PolynomialReduce[x^2, {x - 1}, {x}]  (* univariate division *)
```

### Notes

`PolynomialReduce[poly, {p1, ..., pn}, {x1, ..., xk}]` gives `{{a1, ..., an}, b}`
with `a1 p1 + ... + an pn + b == poly` and `b` fully reduced — no term of the
remainder `b` is divisible by any leading term of the `pi` under the chosen
`MonomialOrder` (default `Lexicographic`). It is the multivariate-division
sibling of `GroebnerBasis` and shares its options (`MonomialOrder`,
`CoefficientDomain`, `Modulus -> p` for `GF(p)`, `ParameterVariables`); free
symbols outside the variable list are treated as coefficient-field parameters. If
the `pi` happen to be a Gröbner basis, `b` is the unique normal form. The
variable list may be omitted, in which case `Variables` is used.
