---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (Kluwer, 1992), ch. 2."
---
### Worked examples

```mathematica
In[1]:= f'[x]
Out[1]= Derivative[1][f][x]
```

```mathematica
In[1]:= D[f[x], x]
Out[1]= Derivative[1][f][x]
```

```mathematica
In[1]:= Derivative[2][Cos]
Out[1]= Derivative[2][Cos]
```

```mathematica
In[1]:= D[f[g[x]], x]
Out[1]= Derivative[1][g][x] Derivative[1][f][g[x]]
```

### Notes

`Derivative[n][f]` is the functional operator representing `f` differentiated `n` times; the surface forms `f'` and `f''` parse to `Derivative[1][f]` and `Derivative[2][f]`. It is the object `D` generates whenever it differentiates an unknown function head, which is why `D[f[x], x]` returns `Derivative[1][f][x]` and the chain rule on `f[g[x]]` yields a product of `Derivative[1]` operators. Note that `Derivative` does not auto-resolve against the known elementary table here: `Derivative[1][Sin]` and `Derivative[2][Cos]` stay in operator form rather than collapsing to `Cos` or `-Cos`. Apply the operator to an explicit argument (e.g. `Derivative[1][f][a]`) to obtain the evaluated-at form.
