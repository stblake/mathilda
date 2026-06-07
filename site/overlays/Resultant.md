---
status: Stable
references:
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\" (3rd ed.), Ch. 6 (resultants and the Sylvester matrix)."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), Ch. 7 (subresultant PRS)."
---
### Worked examples

```mathematica
In[1]:= Resultant[x^2 - 1, x^2 - 4, x]
Out[1]= 9
```

```mathematica
In[1]:= Resultant[x^2 - 2, x^2 - 3, x]
Out[1]= 1
```

```mathematica
In[1]:= Resultant[x^2 + a, x + b, x]
Out[1]= a + b^2
```

```mathematica
In[1]:= Resultant[x^2 - y, x^2 + y, x]
Out[1]= 4 y^2
```

### Notes

`Resultant[p, q, x]` returns a scalar in the remaining variables that
vanishes exactly when `p` and `q` share a common root in `x`. For
`x^2 - 1` and `x^2 - 4` (roots `±1` and `±2`) the value is the nonzero `9`,
confirming no shared root, whereas eliminating `x` from `x^2 - y` and
`x^2 + y` gives `4 y^2`, which is zero only at the shared-root locus `y = 0`.
The computation uses either a Sylvester-matrix determinant or, on the exact
path, a subresultant pseudo-remainder sequence; coefficients may themselves
be polynomials in the other variables, as in `a + b^2`.
