---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), Ch. 3 (monomial extraction from polynomial normal forms)."
---
### Worked examples

```mathematica
In[1]:= Coefficient[x^2 + 3 x + 2, x]
Out[1]= 3
```

```mathematica
In[1]:= Coefficient[x^2 + 3 x + 2, x, 2]
Out[1]= 1
```

```mathematica
In[1]:= Coefficient[a x^2 + b x + c, x, 0]
Out[1]= c
```

```mathematica
In[1]:= Coefficient[3 x^2 y + 2 x y, x, 2]
Out[1]= 3 y
```

### Notes

`Coefficient[expr, form]` returns the coefficient of `form^1` after expanding
`expr`, while the three-argument `Coefficient[expr, form, n]` extracts the
coefficient of `form^n`. Using `n = 0` recovers the term free of `form`, e.g.
the constant `c` in `a x^2 + b x + c`. The extracted coefficient retains any
other variables, so the `x^2` coefficient of `3 x^2 y + 2 x y` is `3 y`. The
`form` is matched structurally against the bases of products, and `n` may be a
non-negative integer (or a rational for Laurent/Puiseux expressions).
