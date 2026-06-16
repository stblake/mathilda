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

```mathematica
In[1]:= Apart[(x^3 + 1)/(x^2 - 1)]
Out[1]= x + 1/(-1 + x)
```

```mathematica
In[1]:= Apart[(2 x + 3)/((x+1)^2 (x^2+1))]
Out[1]= 1/2/(1 + x)^2 + 3/2/(1 + x) + (1 - 3/2 x)/(1 + x^2)
```

```mathematica
In[1]:= Apart[1/(x^2 - 2), Extension -> Sqrt[2]]
Out[1]= -1/2 1/(Sqrt[2] (Sqrt[2] + x)) + 1/2 1/(Sqrt[2] (-Sqrt[2] + x))
```

### Notes

Apart computes the partial-fraction decomposition with respect to the (sole)
variable, the inverse of `Together`. It factors the denominator and splits the
quotient into a sum of terms whose denominators are the factors, including
repeated-factor terms: `1/(x^2 (x+1))` decomposes into `1/x^2 - 1/x + 1/(1+x)`,
recovering both the `1/x^2` and the lower-order `1/x` contributions. Rational
residues appear when the factors are linear over the rationals, as in the
`-1/2` and `3/2` coefficients for `(x+2)/(x^2-1)`. When the numerator degree
meets or exceeds the denominator's, Apart first divides out a polynomial part:
`(x^3 + 1)/(x^2 - 1)` becomes `x + 1/(-1 + x)`. Irreducible quadratic factors
over `Q` are kept intact (the `(1 - 3/2 x)/(1 + x^2)` term), while
`Extension -> Sqrt[2]` splits `x^2 - 2` into the conjugate linear factors
`x +- Sqrt[2]` and produces the corresponding partial fractions over
`Q(Sqrt[2])`.
