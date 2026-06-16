---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (Kluwer, 1992), ch. 2."
---
### Worked examples

```mathematica
In[1]:= Dt[x y]
Out[1]= Dt[x] y + x Dt[y]
```

```mathematica
In[1]:= Dt[Sin[x]]
Out[1]= Cos[x] Dt[x]
```

```mathematica
In[1]:= Dt[Log[x]]
Out[1]= Dt[x]/x
```

```mathematica
In[1]:= Dt[a x, x]
Out[1]= a
```

```mathematica
In[1]:= Dt[x^n]
Out[1]= x^(-1 + n) (n Dt[x] + Dt[n] x Log[x])
```

```mathematica
In[1]:= Dt[f[g[x]]]
Out[1]= Dt[x] Derivative[1][g][x] Derivative[1][f][g[x]]
```

```mathematica
In[1]:= Dt[x^2 y^3]
Out[1]= 2 x Dt[x] y^3 + 3 x^2 y^2 Dt[y]
```

### Notes

`Dt[f]` computes the total differential, treating every symbol as a potential independent variable and emitting `Dt[var]` factors for each one — so `Dt[x y]` gives the full product-rule expansion `Dt[x] y + x Dt[y]`. The two-argument form `Dt[f, x]` is the total derivative with respect to `x`, where other symbols are taken as constants unless they implicitly depend on `x`; `Dt[a x, x]` therefore returns `a`. Elementary functions differentiate through the chain rule with a residual `Dt[x]` factor, as in `Dt[Sin[x]]` and `Dt[Log[x]]`. Constants differentiate to `0` (`Dt[c, x]` gives `0`).
