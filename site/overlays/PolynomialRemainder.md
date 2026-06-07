---
status: Stable
references:
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\" (3rd ed.), Ch. 2 (division with remainder over a field)."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), Ch. 2."
---
### Worked examples

```mathematica
In[1]:= PolynomialRemainder[x^2 - 1, x - 1, x]
Out[1]= 0
```

```mathematica
In[1]:= PolynomialRemainder[x^3 + 2 x^2 + x + 1, x + 1, x]
Out[1]= 1
```

```mathematica
In[1]:= PolynomialRemainder[x^2 + 1, x, x]
Out[1]= 1
```

### Notes

`PolynomialRemainder[p, q, x]` returns the remainder left after dividing `p`
by `q` in `x`; its degree is always strictly less than that of `q`. A zero
remainder, as for `x^2 - 1` divided by `x - 1`, certifies that `q` divides `p`
exactly. Evaluating `p` at a root of a linear divisor recovers the remainder:
`x^3 + 2x^2 + x + 1` at `x = -1` is `1`, matching the output. The companion
`PolynomialQuotient` returns the quotient part; together they satisfy
`p = q*quotient + remainder`. The `Extension -> alpha` option computes the
remainder over `Q(alpha)`.
