---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), on partial fraction decomposition."
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\", on partial fractions and the extended Euclidean algorithm."
---
### Worked examples

```mathematica
In[1]:= Apart[1/(x (x+1))]
Out[1]= 1/x - 1/(1 + x)
```

```mathematica
In[1]:= Apart[(x+2)/(x^2 - 1)]
Out[1]= -1/2/(1 + x) + 3/2/(-1 + x)
```

```mathematica
In[1]:= Apart[1/(x^2 (x+1))]
Out[1]= 1/x^2 - 1/x + 1/(1 + x)
```

### Notes

Apart computes the partial-fraction decomposition with respect to the (sole)
variable, the inverse of `Together`. It factors the denominator and splits the
quotient into a sum of terms whose denominators are the factors, including
repeated-factor terms: `1/(x^2 (x+1))` decomposes into `1/x^2 - 1/x + 1/(1+x)`,
recovering both the `1/x^2` and the lower-order `1/x` contributions. Rational
residues appear when the factors are linear over the rationals, as in the
`-1/2` and `3/2` coefficients for `(x+2)/(x^2-1)`.
