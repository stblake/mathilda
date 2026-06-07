---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), Ch. 3 (multivariate polynomial representation and variable sets)."
---
### Worked examples

```mathematica
In[1]:= Variables[x^2 + y z]
Out[1]= {x, y, z}
```

```mathematica
In[1]:= Variables[a x^2 + b x + c]
Out[1]= {a, b, c, x}
```

```mathematica
In[1]:= Variables[Sin[x] + y]
Out[1]= {Sin[x], y}
```

```mathematica
In[1]:= Variables[x^2 + 3 x + 2]
Out[1]= {x}
```

### Notes

`Variables` collects the independent generators that appear as bases of
non-numeric subexpressions and returns them in canonical sorted order with
duplicates removed. Pure numeric coefficients are ignored, so
`x^2 + 3 x + 2` reports only `{x}`. Compound non-atomic forms are treated as
single generators rather than being broken open: `Sin[x] + y` yields
`{Sin[x], y}`, keeping `Sin[x]` whole. Every symbol that occurs outside
numeric arithmetic is included, which is why parameters like `a, b, c` appear
alongside the main variable `x`.
