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

### Notes

`Dt[f]` computes the total differential, treating every symbol as a potential independent variable and emitting `Dt[var]` factors for each one — so `Dt[x y]` gives the full product-rule expansion `Dt[x] y + x Dt[y]`. The two-argument form `Dt[f, x]` is the total derivative with respect to `x`, where other symbols are taken as constants unless they implicitly depend on `x`; `Dt[a x, x]` therefore returns `a`. Elementary functions differentiate through the chain rule with a residual `Dt[x]` factor, as in `Dt[Sin[x]]` and `Dt[Log[x]]`. Constants differentiate to `0` (`Dt[c, x]` gives `0`).
