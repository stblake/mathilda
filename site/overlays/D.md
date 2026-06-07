---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (Kluwer, 1992), ch. 2."
---
### Worked examples

```mathematica
In[1]:= D[x^n, x]
Out[1]= n x^(-1 + n)
```

```mathematica
In[1]:= D[Exp[a x], x]
Out[1]= a E^(a x)
```

```mathematica
In[1]:= D[x^2 y, {x, 2}]
Out[1]= 2 y
```

```mathematica
In[1]:= D[f[g[x]], x]
Out[1]= Derivative[1][g][x] Derivative[1][f][g[x]]
```

### Notes

`D` is driven by the pattern-based derivative rules in `src/internal/deriv.m`, implementing the chain, product, and quotient rules together with the elementary-function table. Differentiating an unknown function head produces a `Derivative[n][f]` operator rather than evaluating further, so the chain rule on `f[g[x]]` returns a product of such operators. The `{x, n}` form takes the `n`th derivative and treats other symbols as constants by default; use `NonConstants -> {...}` to mark implicit dependencies.
