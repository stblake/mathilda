### Worked examples

```mathematica
In[1]:= Fit[{1, 2, 1.3, 3.75, 2.25}, {1, x}, x]
Out[1]= 0.785 + 0.425 x
```

```mathematica
In[1]:= Fit[{1, 4, 9, 16}, {1, x, x^2}, x]
Out[1]= 0.0 + 0.0 x + 1.0 x^2
```

```mathematica
In[1]:= Fit[{{0, 1}, {1, 2.7}, {2, 7.4}, {3, 20.1}}, {1, x, x^2}, x]
Out[1]= 1.25 - 2.05 x + 2.75 x^2
```

```mathematica
In[1]:= Fit[{{0, 0, 1}, {1, 0, 2}, {0, 1, 3}, {1, 1, 5}}, {1, x, y}, {x, y}]
Out[1]= 0.75 + 1.5 x + 2.5 y
```

### Notes

`Fit[data, {f1, ..., fn}, x]` returns the least-squares linear combination
`a1 f1 + ... + an fn` of the basis functions. Plain `{v1, v2, ...}` data is
taken at abscissae `1, 2, ...`, while `{{x, v}, ...}` pairs supply explicit
abscissae; the perfect-square data recovers `x^2` exactly. The third example
fits a quadratic trend to noisy data, and the last shows a multivariate fit
in two predictors `{x, y}` from `{x, y, v}` rows — the design matrix is
solved by normal equations in either case.
