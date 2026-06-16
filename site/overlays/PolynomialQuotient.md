---
status: Stable
references:
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\" (3rd ed.), Ch. 2 (division with remainder over a field)."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), Ch. 2."
---
### Worked examples

```mathematica
In[1]:= PolynomialQuotient[x^2 - 1, x - 1, x]
Out[1]= 1 + x
```

```mathematica
In[1]:= PolynomialQuotient[x^3 + 2 x^2 + x + 1, x + 1, x]
Out[1]= x + x^2
```

```mathematica
In[1]:= PolynomialQuotient[x^2 + 1, x, x]
Out[1]= x
```

```mathematica
In[1]:= PolynomialQuotient[x^4 - 2, x^2 - Sqrt[2], x, Extension -> Sqrt[2]]
Out[1]= Sqrt[2] + x^2
```

### Notes

`PolynomialQuotient[p, q, x]` performs long division of `p` by `q` in the
variable `x` and returns only the quotient, discarding the remainder. So
`(x^2 + 1) / x` yields `x` (dropping the `1` remainder) and `x^2 - 1` divides
exactly by `x - 1` to give `1 + x`. Pair it with `PolynomialRemainder` to
recover the full relation `p = q*quotient + remainder`. The `Extension -> alpha`
option divides over `Q(alpha)`; the default `None`/`Automatic` divides over
the rationals.
